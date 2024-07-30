from latch_asgi.context.http import Context, HandlerResult
from latch_asgi.framework.http import HTTPInternalServerError
from latch_o11y.o11y import trace_function_with_span
from opentelemetry.trace import Span, get_tracer

from ..entrypoint import k_proc, start_kernel_proc

tracer = get_tracer(__name__)


@trace_function_with_span(tracer)
async def readyz(s: Span, ctx: Context) -> HandlerResult:
    proc = k_proc.proc
    if proc is None:
        raise HTTPInternalServerError("Kernel process not started")

    # note: if the kernel proc died for some reason
    if proc.returncode is not None:
        print("Restarting kernel process")
        await start_kernel_proc()

    return "OK"
