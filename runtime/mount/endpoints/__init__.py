from latch_asgi.context.http import Route
from latch_asgi.context.websocket import Route as WsRoute

from .readyz import readyz
from .run import run

http_routes: dict[str, Route] = {"/readyz": (["GET"], readyz)}

websocket_routes: dict[str, WsRoute] = {"/run": run}
