import asyncio
import concurrent.futures as cf
import threading
from collections.abc import Coroutine
from socket import socket
from threading import Thread
from typing import Any, TypeVar

from socketio import SocketIo

T = TypeVar("T")


class SocketIoThread(Thread):
    loop: asyncio.AbstractEventLoop | None = None
    conn: SocketIo | None = None
    lock: asyncio.Lock | None = None

    def __init__(self, *, socket: socket) -> None:
        super().__init__(name="socket_io")

        self.socket = socket
        self.shutdown = asyncio.Event()
        self.initialized = threading.Event()
        self.lock = asyncio.Lock()

    def run(self) -> None:
        async def f():
            self.loop = asyncio.get_running_loop()
            self.conn = await SocketIo.from_socket(self.socket)

            self.initialized.set()
            await self.shutdown.wait()

        asyncio.run(f())

    def call_fut(self, coro: Coroutine[None, None, T]) -> cf.Future[T]:
        assert self.loop is not None

        return asyncio.run_coroutine_threadsafe(coro, self.loop)

    def send_fut(self, data: object) -> cf.Future[None]:
        assert self.conn is not None

        return self.call_fut(self.conn.send(data))

    async def send(self, data: object) -> None:
        async with self.lock:
            await asyncio.wrap_future(self.send_fut(data))

    def recv_fut(self) -> cf.Future[Any]:
        assert self.conn is not None

        return self.call_fut(self.conn.recv())

    async def recv(self) -> Any:
        return await asyncio.wrap_future(self.recv_fut())
