import asyncio
from dataclasses import dataclass
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
        self.conn.send_fut({
            "type": "kernel_stdio",
            "active_cell": self.kernel.active_cell,
            "stream": self.name,
            # todo(maximsmol): this is a bit silly because we are going to have
            # a TextIOWrapper above that just encoded this for us
            "data": data,
        }).result()
        return len(b)


def text_socket_writer(x: SocketWriter) -> TextIOWrapper:
    return TextIOWrapper(
        BufferedWriter(x), encoding="utf-8", errors="replace", line_buffering=True
    )
