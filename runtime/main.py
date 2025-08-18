import sys

from latch_o11y.o11y import setup as setup_o11y

setup_o11y()

import asyncio
import signal

import uvicorn
from latch_asgi.server import LatchASGIServer

from .config import config
from .mount.endpoints import http_routes, websocket_routes
from .mount.entrypoint import shutdown, start_kernel_proc

latch_server = LatchASGIServer(
    http_routes=http_routes,
    websocket_routes=websocket_routes,
    startup_tasks=[start_kernel_proc()],
    shutdown_tasks=[shutdown()],
)


async def app(scope, receive, send):
    await latch_server.raw_app(scope, receive, send)


async def _run_uvicorn() -> None:
    uv_config = uvicorn.Config(
        app,
        host="0.0.0.0",
        port=5000,
        ws="websockets",
        reload=config.auto_reload,
        log_config=None,
    )
    server = uvicorn.Server(uv_config)

    shutdown_event = asyncio.Event()

    def shutdown_signal(*_args: object) -> None:
        shutdown_event.set()

    loop = asyncio.get_running_loop()
    loop.add_signal_handler(signal.SIGTERM, shutdown_signal)
    loop.add_signal_handler(signal.SIGINT, shutdown_signal)

    async def _shutdown_watcher() -> None:
        await shutdown_event.wait()
        server.should_exit = True

    await asyncio.gather(server.serve(), _shutdown_watcher())


if __name__ == "__main__":
    if sys.platform == "linux":
        from ctypes import CDLL

        libc = CDLL("libc.so.6")
        PR_SET_NAME = 15  # https://github.com/torvalds/linux/blob/2df0c02dab829dd89360d98a8a1abaa026ef5798/include/uapi/linux/prctl.h#L56
        libc.prctl(PR_SET_NAME, b"latch-plots")

    import uvloop

    uvloop.install()
    asyncio.run(_run_uvicorn())
