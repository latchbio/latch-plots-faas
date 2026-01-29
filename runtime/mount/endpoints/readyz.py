from latch_asgi.context.http import Context, HandlerResult
from latch_asgi.framework.http import HTTPInternalServerError

from ..entrypoint import main_k_proc, start_kernel_proc


async def readyz(ctx: Context) -> HandlerResult:
    proc = main_k_proc.proc
    if proc is None:
        raise HTTPInternalServerError("Kernel process not started")

    # note: if the kernel proc died for some reason
    if proc.returncode is not None:
        print("Restarting kernel process")
        await start_kernel_proc(main_k_proc)

    return "OK"
