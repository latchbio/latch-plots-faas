import asyncio
import concurrent.futures as cf
import contextlib
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
    send_queue: asyncio.Queue[tuple[object, cf.Future[None]]] | None = None
    sender_task: asyncio.Task[None] | None = None

    def __init__(self, *, socket: socket) -> None:
        super().__init__(name="socket_io")

        self.socket = socket
        self.shutdown = asyncio.Event()
        self.initialized = threading.Event()

    async def _send_loop(self) -> None:
        assert self.conn is not None
        assert self.send_queue is not None

        while True:
            data, future = await self.send_queue.get()
            try:
                await self.conn.send(data)
                future.set_result(None)
            except Exception as e:
                future.set_exception(e)
            finally:
                self.send_queue.task_done()

    def run(self) -> None:
        async def f() -> None:
            self.loop = asyncio.get_running_loop()
            self.conn = await SocketIo.from_socket(self.socket)
            self.send_queue = asyncio.Queue()
            self.sender_task = asyncio.create_task(self._send_loop())

            self.initialized.set()
            await self.shutdown.wait()

            if self.sender_task is not None:
                self.sender_task.cancel()
                with contextlib.suppress(asyncio.CancelledError):
                    await self.sender_task

        asyncio.run(f())

    def call_fut(self, coro: Coroutine[None, None, T]) -> cf.Future[T]:
        assert self.loop is not None

        return asyncio.run_coroutine_threadsafe(coro, self.loop)

    def send_fut(self, data: object) -> cf.Future[None]:
        assert self.loop is not None
        assert self.send_queue is not None
        future: cf.Future[None] = cf.Future()

        self.send_queue.put_nowait((data, future))

        return future

    async def send(self, data: object) -> None:
        await asyncio.wrap_future(self.send_fut(data))

    def recv_fut(self) -> cf.Future[Any]:
        assert self.conn is not None

        return self.call_fut(self.conn.recv())

    async def recv(self) -> Any:
        return await asyncio.wrap_future(self.recv_fut())
