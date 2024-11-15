import asyncio
import concurrent.futures as cf
import threading
from collections.abc import Coroutine
from threading import Thread
from typing import TypeVar

from socketio import SocketIo

T = TypeVar("T")


class WidgetUpdateThread(Thread):
    loop: asyncio.AbstractEventLoop | None = None
    conn: SocketIo | None = None

    def __init__(self) -> None:
        super().__init__(name="widget_update")

        self.shutdown = asyncio.Event()
        self.initialized = threading.Event()

    def run(self) -> None:
        async def f():
            self.loop = asyncio.get_running_loop()
            self.initialized.set()
            await self.shutdown.wait()

        asyncio.run(f())

    def send(self, coro: Coroutine[None, None, T]) -> cf.Future[T]:
        assert self.loop is not None

        return asyncio.run_coroutine_threadsafe(coro, self.loop)
