import asyncio
import threading
from dataclasses import dataclass, field
from io import BufferedWriter, RawIOBase, TextIOWrapper, UnsupportedOperation
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


@dataclass(kw_only=True)
class SocketWriter(RawIOBase):
    conn: SocketIoThread
    kernel: "Kernel"
    name: str
    loop: asyncio.AbstractEventLoop
    buffer: list[tuple[int | None, str, str]] = field(default_factory=list, init=False)
    buffer_lock: threading.Lock = field(default_factory=threading.Lock, init=False)
    flush_task: asyncio.Task | None = field(default=None, init=False)

    def __post_init__(self) -> None:
        self.buffer = []
        self.buffer_lock = threading.Lock()
        self.flush_task = asyncio.run_coroutine_threadsafe(
            self._periodic_flush(), self.loop
        )

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

        with self.buffer_lock:
            self.buffer.append((self.kernel.active_cell, self.name, data))

        return len(b)

    def _flush_buffer(self) -> None:
        """Flush all buffered messages to the socket."""
        with self.buffer_lock:
            messages = self.buffer.copy()
            self.buffer.clear()

        for active_cell, stream, data in messages:
            self.conn.send_fut({
                "type": "kernel_stdio",
                "active_cell": active_cell,
                "stream": stream,
                "data": data,
            }).result()

    async def _periodic_flush(self) -> None:
        """Periodically flush the buffer every 0.5 seconds."""
        try:
            while True:
                await asyncio.sleep(0.5)
                self._flush_buffer()
        except asyncio.CancelledError:
            # Task was cancelled, exit gracefully
            pass

    @override
    def close(self) -> None:
        """Close the writer and flush any remaining buffered data."""
        if self.flush_task is not None:
            self.flush_task.cancel()
            try:
                self.flush_task.result(timeout=1.0)
            except Exception:
                pass

        # Final flush to ensure no data is lost
        self._flush_buffer()

        super().close()


def text_socket_writer(x: SocketWriter) -> TextIOWrapper:
    return TextIOWrapper(
        BufferedWriter(x), encoding="utf-8", errors="replace", line_buffering=True
    )
