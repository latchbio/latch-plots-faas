from latch_asgi.context.http import Context, HandlerResult
from latch_asgi.framework.http import HTTPInternalServerError

from ..entrypoint import k_proc


async def readyz(ctx: Context) -> HandlerResult:
    proc = k_proc.proc
    if proc is None:
        raise HTTPInternalServerError("Kernel process not started")

    return "OK"
