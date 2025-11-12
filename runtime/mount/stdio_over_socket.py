import asyncio
import concurrent.futures as cf
import threading
from dataclasses import dataclass, field
from io import BufferedWriter, RawIOBase, TextIOWrapper, UnsupportedOperation
from itertools import groupby
from operator import itemgetter
from typing import TYPE_CHECKING

from socketio_thread import SocketIoThread
from typing_extensions import override

if TYPE_CHECKING:
    from _typeshed import ReadableBuffer
    from kernel import Kernel

# todo(kenny): without snapshot these warnings occurred before connection was
# established. Snapshot load pushes warnings to after connection and cause
# console to show errors. Actual fix is addressing root cause in each
# dependency but this might break code for customers.
warning_msgs = [
        "/opt/mamba/envs/plots-faas/lib/python3.11/site-packages/dask/dataframe/__init__.py:31: FutureWarning:",
        "/opt/mamba/envs/plots-faas/lib/python3.11/site-packages/numba/core/decorators.py:246: RuntimeWarning:",
        "/opt/mamba/envs/plots-faas/lib/python3.11/site-packages/anndata/utils.py:429: FutureWarning:",
        ]


flush_interval = 0.5


@dataclass(kw_only=True)
class SocketWriter(RawIOBase):
    conn: SocketIoThread
    kernel: "Kernel"
    name: str
    loop: asyncio.AbstractEventLoop
    _buffer: list[tuple[str, str | None]] = field(default_factory=list, init=False)
    _buffer_lock: threading.Lock = field(default_factory=threading.Lock, init=False)
    _flusher_task: cf.Future[None] | None = field(default=None, init=False)

    def _start_flusher(self) -> None:
        if self._flusher_task is None:
            self._flusher_task = asyncio.run_coroutine_threadsafe(
                self._flush_loop(), self.conn.loop
            )

    async def _flush_loop(self) -> None:
        while True:
            await asyncio.sleep(flush_interval)
            await self._flush()

    async def _flush(self) -> None:
        with self._buffer_lock:
            if len(self._buffer) == 0:
                return

            items = list(self._buffer)
            self._buffer.clear()

        for cell_id, group in groupby(items, key=itemgetter(1)):
            combined_data = "".join(data for data, _ in group)
            await self.conn.send({
                "type": "kernel_stdio",
                "active_cell": cell_id,
                "stream": self.name,
                "data": combined_data,
            })

    @override
    def fileno(self) -> int:
        # todo(maximsmol): if this causes problems, allow a workaround of some sort here
        # either return the default fd for stdout/stderr if it's not actually important
        # or we need to make a temporary file and read from it as a fallback for things
        # that want to reopen the stream for some reason
        raise UnsupportedOperation(
            "sys.stdout.fileno() is not supported in notebooks: "
            "data is redirected over the kernel message stream "
            "and does not have a corresponding file descriptor"
        )

    @override
    def writable(self) -> bool:
        return True

    @override
    def write(self, __b: "ReadableBuffer") -> int | None:
        b = bytes(__b)
        data = b.decode(errors="replace")
        for x in warning_msgs:
            if x in data:
                return len(b)

        if self._flusher_task is None:
            self._start_flusher()

        with self._buffer_lock:
            self._buffer.append((data, self.kernel.active_cell))

        return len(b)


def text_socket_writer(x: SocketWriter) -> TextIOWrapper:
    return TextIOWrapper(
        BufferedWriter(x), encoding="utf-8", errors="replace", line_buffering=True
    )
