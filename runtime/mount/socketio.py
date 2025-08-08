import asyncio
import socket
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Self

import orjson

sys.path.append(str(Path(__file__).parent.absolute()))
from utils import orjson_encoder

logfile = Path("/tmp/socketio_messages.log")  # noqa: S108
logfile.touch()


def append_log_message(message: str) -> None:
    with logfile.open("a") as f:
        f.write(message)


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
        append_log_message(f"sending message of {len(data)} bytes")

        async with self.wlock:
            await self.loop.sock_sendall(self.sock, header)
            await self.loop.sock_sendall(self.sock, data)

    async def send(self, data: object) -> None:
        await self.send_bytes(
            orjson.dumps(
                data, option=orjson.OPT_SERIALIZE_NUMPY, default=orjson_encoder
            )
        )

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

        return orjson.loads(data)
