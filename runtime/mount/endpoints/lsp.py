import asyncio
import json
from typing import Any

from latch_asgi.context.websocket import Context, HandlerResult
from latch_asgi.framework.websocket import WebsocketConnectionClosedError, receive_json
from latch_o11y.o11y import trace_app_function_with_span
from opentelemetry.trace import Span


async def send_lsp_msg(stdin: asyncio.StreamWriter, message: dict[str, Any]) -> None:
    global request_id
    request_id += 1

    body = json.dumps(message, indent=2)
    content_length = len(body)
    header = f"Content-Length: {content_length}\r\n\r\n"

    stdin.write((header + body).encode("utf-8"))
    await stdin.drain()


async def recv_lsp_msg(stdout: asyncio.StreamReader) -> dict[str, Any]:
    response = await stdout.readuntil(b"\r\n\r\n")

    content_length = 0
    headers = response.decode("utf-8").split("\r\n")
    for header in headers:
        if header.startswith("Content-Length:"):
            content_length = int(header.split(": ")[1])
            break

    body = await stdout.readexactly(content_length)
    return json.loads(body.decode("utf-8"))


@trace_app_function_with_span
async def lsp_proxy(s: Span, ctx: Context) -> HandlerResult:
    # todo(rteqs): auth

    async def poll_lsp_msg(stdout: asyncio.StreamReader) -> None:
        while True:
            # todo(rteqs): probably need to include header as well
            data = await recv_lsp_msg(stdout)
            await ctx.send_message(data)

    async def poll_lsp_err(stderr: asyncio.StreamReader) -> None:
        while True:
            line = (await stderr.readline()).decode("utf-8")
            await ctx.send_message(line)

    stdin_pipe = asyncio.subprocess.PIPE
    stdout_pipe = asyncio.subprocess.PIPE
    stderr_pipe = asyncio.subprocess.PIPE

    proc = await asyncio.subprocess.create_subprocess_exec(
        "/opt/mamba/envs/plots-faas/bin/pyright-langserver",
        "--stdio",
        stdin=stdin_pipe,
        stdout=stdout_pipe,
        stderr=stderr_pipe,
    )
    stdin = proc.stdin
    stdout = proc.stdout
    stderr = proc.stderr

    try:
        async with asyncio.TaskGroup() as tg:
            tg.create_task(poll_lsp_msg(stdout))
            tg.create_task(poll_lsp_err(stderr))
            while True:
                msg = await receive_json(ctx.receive)
                await send_lsp_msg(stdin, msg)

    except WebsocketConnectionClosedError:
        ...

    proc.kill()
    await proc.wait()
    return
