import asyncio
import json

from latch_asgi.context.websocket import Context, HandlerResult
from latch_asgi.framework.websocket import WebsocketConnectionClosedError, receive_json
from latch_o11y.o11y import trace_app_function_with_span
from opentelemetry.trace import Span


async def recv_lsp_msg(stdout: asyncio.StreamReader) -> bytes:
    headers_raw = await stdout.readuntil(b"\r\n\r\n")

    content_length = 0
    headers_parsed = headers_raw.decode("utf-8").split("\r\n")
    for header in headers_parsed:
        if header.startswith("Content-Length:"):
            content_length = int(header.split(": ")[1])
            break

    body = await stdout.readexactly(content_length)
    return headers_raw + body


@trace_app_function_with_span
async def lsp_proxy(s: Span, ctx: Context) -> HandlerResult:
    # todo(rteqs): auth
    await ctx.accept_connection()

    async def poll_lsp_msg(stdout: asyncio.StreamReader) -> None:
        while True:
            # todo(rteqs): probably need to include header as well
            data = await recv_lsp_msg(stdout)
            await ctx.send_message(json.dumps({"type": "msg", "data": data}))

    async def poll_lsp_err(stderr: asyncio.StreamReader) -> None:
        while True:
            line = (await stderr.readline()).decode("utf-8")
            # todo(rteqs): probably don't need to send cz server should handle everything including restarting pyright if it fails
            await ctx.send_message(json.dumps({"type": "error", "data": line}))

    stdin_pipe = asyncio.subprocess.PIPE
    stdout_pipe = asyncio.subprocess.PIPE
    stderr_pipe = asyncio.subprocess.PIPE

    proc = await asyncio.subprocess.create_subprocess_exec(
        "/opt/mamba/envs/plots-faas/bin/node",
        "/opt/mamba/envs/plots-faas/bin/pyright-langserver",
        "--stdio",
        stdin=stdin_pipe,
        stdout=stdout_pipe,
        stderr=stderr_pipe,
    )
    stdin = proc.stdin
    stdout = proc.stdout
    stderr = proc.stderr

    assert stdin is not None
    assert stdout is not None
    assert stderr is not None

    try:
        async with asyncio.TaskGroup() as tg:
            tg.create_task(poll_lsp_msg(stdout))
            tg.create_task(poll_lsp_err(stderr))
            while True:
                # todo(rteqs): probably need different types for this
                msg = await receive_json(ctx.receive)
                stdin.write(msg)
                await stdin.drain()

    except WebsocketConnectionClosedError:
        ...

    proc.kill()
    await proc.wait()
    return "Ok"
