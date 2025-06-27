import asyncio
from dataclasses import dataclass
from io import BufferedWriter, RawIOBase, TextIOWrapper, UnsupportedOperation
from typing import TYPE_CHECKING

from socketio_thread import SocketIoThread
from typing_extensions import override

if TYPE_CHECKING:
    from _typeshed import ReadableBuffer
    from kernel import Kernel


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
        self.conn.send_fut({
            "type": "kernel_stdio",
            "active_cell": self.kernel.active_cell,
            # "snapshot_status": self.kernel.snapshot_status,
            "stream": self.name,
            # todo(maximsmol): this is a bit silly because we are going to have
            # a TextIOWrapper above that just encoded this for us
            "data": b.decode(errors="replace"),
        }).result()
        return len(b)


def text_socket_writer(x: SocketWriter) -> TextIOWrapper:
    return TextIOWrapper(
        BufferedWriter(x), encoding="utf-8", errors="replace", line_buffering=True
    )
