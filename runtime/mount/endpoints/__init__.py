from latch_asgi.context.http import Route
from latch_asgi.context.websocket import Route as WsRoute

from .agent import agent
from .headless import spawn_headless_browser, stop_headless_browser_endpoint
from .lsp import lsp_proxy
from .readyz import readyz
from .run import run

http_routes: dict[str, Route] = {
    "/readyz": (["GET"], readyz),
    "/headless/spawn": (["POST"], spawn_headless_browser),
    "/headless/stop": (["POST"], stop_headless_browser_endpoint),
}

websocket_routes: dict[str, WsRoute] = {"/run": run, "/lsp": lsp_proxy, "/agent": agent}
