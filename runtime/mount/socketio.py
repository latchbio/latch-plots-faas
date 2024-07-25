import asyncio
import base64
import hashlib
import json
import re
import socket
import struct
from dataclasses import dataclass
from enum import Enum
from typing import Any, Self


@dataclass(kw_only=True)
class SocketIo:
    sock: socket.socket
    loop: asyncio.AbstractEventLoop

    rlock: asyncio.Lock
    wlock: asyncio.Lock

    @classmethod
    async def from_socket(cls, sock: socket.socket) -> Self:
        loop = asyncio.get_event_loop()

        return cls(sock=sock, loop=loop, rlock=asyncio.Lock(), wlock=asyncio.Lock())

    async def send_bytes(self, data: bytes) -> None:
        header = struct.pack("<q", len(data))

        async with self.wlock:
            await self.loop.sock_sendall(self.sock, header)
            await self.loop.sock_sendall(self.sock, data)

    async def send(self, data: object) -> None:
        await self.send_bytes(json.dumps(data, allow_nan=False).encode())

    async def recv(self) -> Any:
        async with self.rlock:
            header = await self.loop.sock_recv(self.sock, 8)
            if len(header) == 0:
                raise EOFError

            (l,) = struct.unpack("<q", header)
            assert isinstance(l, int)

            data = b""
            while l > 0:
                chunk = await self.loop.sock_recv(self.sock, l)
                l -= len(chunk)
                data += chunk

        return json.loads(data)


class Opcode(Enum):
    continuation = 0x0
    text = 0x1
    binary = 0x2
    close = 0x8


class WebsocketClosedError(RuntimeError):
    pass


# https://developer.mozilla.org/en-US/docs/Web/API/WebSockets_API/Writing_WebSocket_servers
# todo(rteqs): this only has the stuff we need, not the full spec
@dataclass(kw_only=True)
class WebSocketIo:
    sock: socket.socket
    loop: asyncio.AbstractEventLoop

    rlock: asyncio.Lock
    wlock: asyncio.Lock

    @classmethod
    async def from_socket(cls, sock: socket.socket) -> Self:
        loop = asyncio.get_event_loop()

        return cls(sock=sock, loop=loop, rlock=asyncio.Lock(), wlock=asyncio.Lock())

    async def send_bytes(self, data: bytes, opcode: Opcode = Opcode.binary) -> None:
        payload_length = len(data)

        frame_header = bytearray()
        frame_header.append(0x80 | opcode.value)

        if payload_length <= 125:
            frame_header.append(payload_length)
        elif payload_length <= 65535:
            frame_header.append(126)
            frame_header.extend(struct.pack(">H", payload_length))
        else:
            frame_header.append(127)
            frame_header.extend(struct.pack(">Q", payload_length))

        payload = frame_header + data
        async with self.wlock:
            await self.loop.sock_sendall(self.sock, payload)

    async def send(self, data: object) -> None:
        await self.send_bytes(json.dumps(data, allow_nan=False).encode())

    async def recv(self) -> Any:
        data = b""
        fin = False

        while not fin:
            frame_header = await self.loop.sock_recv(self.sock, 2)
            if len(frame_header) == 0:
                break

            fin_and_opcode, mask_and_payload_len = frame_header
            fin = (fin_and_opcode & 0x80) != 0
            opcode = fin_and_opcode & 0xF

            if opcode == Opcode.close.value:
                raise WebsocketClosedError

            payload_length = mask_and_payload_len & 0x7F
            if payload_length == 126:
                extended_payload_length = await self.loop.sock_recv(self.sock, 2)
                payload_length = struct.unpack(">H", extended_payload_length)[0]
            elif payload_length == 127:
                extended_payload_length = await self.loop.sock_recv(self.sock, 8)
                payload_length = struct.unpack(">Q", extended_payload_length)[0]

            mask = (mask_and_payload_len & 0x80) != 0
            masking_key = await self.loop.sock_recv(self.sock, 4) if mask else None

            payload_data = await self.loop.sock_recv(self.sock, payload_length)

            if mask and masking_key is not None:
                payload = bytearray(payload_data)
                for i in range(len(payload)):
                    payload[i] ^= masking_key[i % 4]
                payload_data = bytes(payload)

            data += payload_data

        return json.loads(data)

    async def close(self, code: int = 1000) -> None:
        close_code = struct.pack(">H", code)
        await self.send_bytes(data=close_code, opcode=Opcode.close)


upgrade_pattern = re.compile(rb"Upgrade:\s+websocket")
sec_websocket_key_pattern = re.compile(
    rb"Sec-WebSocket-Key:\s+(?P<key>[+/0-9A-Za-z=]+)"
)


async def upgrade_connection(sock: socket.socket) -> None:
    loop = asyncio.get_running_loop()
    data = await loop.sock_recv(sock, 4096)

    upgrade_match = upgrade_pattern.search(data)
    sec_websocket_key_match = sec_websocket_key_pattern.search(data)

    if not upgrade_match or not sec_websocket_key_match:
        raise ValueError("Invalid WebSocket upgrade request")

    websocket_key = sec_websocket_key_match.group("key").decode("utf-8")

    magic_string = "258EAFA5-E914-47DA-95CA-C5AB0DC85B11"
    accept_key = base64.b64encode(
        hashlib.sha1((websocket_key + magic_string).encode("utf-8")).digest()  # noqa: S324
    ).decode("utf-8")

    res = (
        "HTTP/1.1 101 Switching Protocols\r\n"
        "Upgrade: websocket\r\n"
        "Connection: Upgrade\r\n"
        f"Sec-WebSocket-Accept: {accept_key}\r\n\r\n"
    )

    await loop.sock_sendall(sock, res.encode("utf-8"))
