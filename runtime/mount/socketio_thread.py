import asyncio
from socket import socket
from threading import Thread
import concurrent.futures as cf
import threading
from typing import Any, TypeVar
from collections.abc import Coroutine

from socketio import SocketIo

T = TypeVar("T")


class SocketIoThread(Thread):
    loop: asyncio.AbstractEventLoop | None = None
    conn: SocketIo | None = None

    def __init__(self, *, socket: socket) -> None:
        super().__init__()

        self.socket = socket
        self.shutdown = asyncio.Event()
        self.initialized = threading.Event()

    def run(self) -> None:
        async def f():
            self.loop = asyncio.get_running_loop()
            self.conn = await SocketIo.from_socket(self.socket)

            self.initialized.set()
            await self.shutdown.wait()

        asyncio.run(f())

    def _call(self, coro: Coroutine[None, None, T]) -> cf.Future[T]:
        assert self.loop is not None

        return asyncio.run_coroutine_threadsafe(coro, self.loop)

    def send_fut(self, data: object) -> cf.Future[None]:
        assert self.conn is not None

        return self._call(self.conn.send(data))

    async def send(self, data: object) -> None:
        await asyncio.wrap_future(self.send_fut(data))

    def recv_fut(self) -> cf.Future[Any]:
        assert self.conn is not None

        return self._call(self.conn.recv())

    async def recv(self) -> Any:
        return await asyncio.wrap_future(self.recv_fut())
