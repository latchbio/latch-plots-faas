import asyncio
import gzip
import socket
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Self

import orjson

sys.path.append(str(Path(__file__).parent.absolute()))
from utils import orjson_encoder


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
        data = gzip.compress(data, compresslevel=0)
        header = struct.pack("<q", len(data))

        async with self.wlock:
            await self.loop.sock_sendall(self.sock, header)
            await self.loop.sock_sendall(self.sock, data)

    async def send(self, data: object) -> None:
        await self.send_bytes(
            orjson.dumps(
                data, option=orjson.OPT_SERIALIZE_NUMPY, default=orjson_encoder
            )
        )

    async def recv_raw(self) -> bytes:
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

        return gzip.decompress(data)

    async def recv(self) -> Any:
        return orjson.loads(await self.recv_raw())
