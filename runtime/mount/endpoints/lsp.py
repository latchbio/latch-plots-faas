import asyncio
import os

import orjson
from latch_asgi.context.websocket import Context, HandlerResult
from latch_asgi.framework.websocket import WebsocketConnectionClosedError, receive_data
from latch_o11y.o11y import trace_app_function_with_span
from opentelemetry.trace import Span


async def send_lsp_msg(stdin: asyncio.StreamWriter, msg: bytes) -> None:
    content_length = len(msg)
    header = f"Content-Length: {content_length}\r\n\r\n".encode("ascii")
    stdin.write(header + msg)
    await stdin.drain()


async def recv_lsp_msg(stdout: asyncio.StreamReader) -> bytes:
    headers_raw = await stdout.readuntil(b"\r\n\r\n")

    content_length = 0
    headers_parsed = headers_raw.decode("ascii").split("\r\n")
    for header in headers_parsed:
        if header.startswith("Content-Length:"):
            content_length = int(header.split(": ")[1])
            break

    return await stdout.readexactly(content_length)


@trace_app_function_with_span
async def lsp_proxy(s: Span, ctx: Context) -> HandlerResult:
    await ctx.accept_connection()

    async def poll_lsp_msg(stdout: asyncio.StreamReader) -> None:
        while True:
            data = await recv_lsp_msg(stdout)
            await ctx.send_message(data.decode("utf-8"))

    async def poll_lsp_err(stderr: asyncio.StreamReader) -> None:
        while True:
            line = (await stderr.readline()).decode("utf-8")
            print(line)
            # todo(rteqs): error handling

    stdin_pipe = asyncio.subprocess.PIPE
    stdout_pipe = asyncio.subprocess.PIPE
    stderr_pipe = asyncio.subprocess.PIPE

    proc = await asyncio.subprocess.create_subprocess_exec(
        "/opt/mamba/envs/plots-faas/bin/pyright-langserver",
        "--stdio",
        stdin=stdin_pipe,
        stdout=stdout_pipe,
        stderr=stderr_pipe,
        env={"PATH": "/opt/mamba/envs/plots-faas/bin:" + os.environ["PATH"]},
    )
    stdin = proc.stdin
    stdout = proc.stdout
    stderr = proc.stderr

    assert stdin is not None
    assert stdout is not None
    assert stderr is not None

    async with asyncio.TaskGroup() as tg:
        poll_msg_task = tg.create_task(poll_lsp_msg(stdout))
        poll_err_task = tg.create_task(poll_lsp_err(stderr))
        try:
            while True:
                msg = await receive_data(ctx.receive)
                if isinstance(msg, str):
                    msg = msg.encode("utf-8")
                await send_lsp_msg(stdin, msg)

        except WebsocketConnectionClosedError:
            ...

        finally:
            poll_msg_task.cancel()
            poll_err_task.cancel()

    shutdown_msg = orjson.dumps({"jsonrpc": "2.0", "method": "shutdown"})
    await send_lsp_msg(stdin, shutdown_msg)

    exit_msg = orjson.dumps({"jsonrpc": "2.0", "method": "exit"})
    await send_lsp_msg(stdin, exit_msg)

    proc.kill()
    await proc.wait()
    return "Ok"
