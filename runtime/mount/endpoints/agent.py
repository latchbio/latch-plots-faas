import asyncio

import orjson
from latch_asgi.context.websocket import Context, HandlerResult
from latch_asgi.framework.websocket import (
    receive_json,
)
from latch_o11y.o11y import trace_app_function_with_span
from opentelemetry.trace import Span

from ..entrypoint import a_proc, pod_id, pod_session_id, start_agent_proc

connection_idx = 0
agent_start_lock = asyncio.Lock()


@trace_app_function_with_span
async def agent(s: Span, ctx: Context) -> HandlerResult:
    global connection_idx

    await ctx.accept_connection()

    s.set_attributes({
        "pod_id": pod_id,
        "pod_session_id": pod_session_id,
        "connection_idx": connection_idx,
    })

    async with agent_start_lock:
        if a_proc.conn_a is None:
            await start_agent_proc()

    conn_a = a_proc.conn_a
    if conn_a is None:
        await ctx.send_message(
            orjson.dumps({
                "type": "agent_error",
                "error": "Agent process failed to start",
                "fatal": True
            }).decode()
        )
        return "Agent not available"

    await conn_a.send({
        "type": "init",
    })

    await ctx.send_message(
        orjson.dumps({
            "type": "agent_ready",
            "connection_idx": connection_idx,
        }).decode()
    )

    connection_idx += 1

    async def forward_to_agent() -> dict:
        while True:
            msg = await receive_json(ctx.receive)
            await conn_a.send(msg)

    async def forward_to_frontend() -> dict:
        while True:
            msg = await conn_a.recv()
            await ctx.send_message(orjson.dumps(msg).decode())

    await asyncio.gather(
        forward_to_agent(),
        forward_to_frontend(),
        return_exceptions=True
    )

    return "Ok"
