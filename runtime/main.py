from latch_o11y.o11y import setup as setup_o11y

setup_o11y()

import asyncio

from hypercorn.asyncio import serve
from hypercorn.config import Config as HypercornConfig
from latch_asgi.server import LatchASGIServer

from .config import config
from .mount.endpoints import http_routes, websocket_routes
from .mount.entrypoint import shutdown, start_kernel_proc

cfg = HypercornConfig()
cfg.bind = ["0.0.0.0:5000"]
cfg.use_reloader = config.auto_reload
cfg.graceful_timeout = 0.1


latch_server = LatchASGIServer(
    http_routes=http_routes,
    websocket_routes=websocket_routes,
    startup_tasks=[start_kernel_proc()],
    shutdown_tasks=[shutdown()],
)

app = latch_server.raw_app

if __name__ == "__main__":
    import signal

    import uvloop

    uvloop.install()
    loop = asyncio.new_event_loop()

    shutdown_event = asyncio.Event()

    async def await_shutdown() -> None:
        await shutdown_event.wait()

    def shutdown_signal(*args) -> None:
        shutdown_event.set()

    loop.add_signal_handler(signal.SIGTERM, shutdown_signal)
    loop.add_signal_handler(signal.SIGINT, shutdown_signal)
    loop.run_until_complete(serve(app, cfg, shutdown_trigger=await_shutdown))
