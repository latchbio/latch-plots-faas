import asyncio

import orjson
from latch_asgi.context.websocket import Context, HandlerResult
from latch_asgi.framework.websocket import (
    receive_json,
)
from latch_o11y.o11y import trace_app_function_with_span
from opentelemetry.trace import Span

import runtime.mount.entrypoint as entrypoint_module

from ..entrypoint import a_proc, pod_id, pod_session_id, start_agent_proc

connection_idx = 0
agent_start_lock = asyncio.Lock()


@trace_app_function_with_span
async def agent(s: Span, ctx: Context) -> HandlerResult:
    global connection_idx

    print(f"[agent-endpoint] New agent WebSocket connection (idx={connection_idx})")
    await ctx.accept_connection()
    print(f"[agent-endpoint] Connection accepted (idx={connection_idx})")

    s.set_attributes({
        "pod_id": pod_id,
        "pod_session_id": pod_session_id,
        "connection_idx": connection_idx,
    })

    entrypoint_module.current_agent_ctx = ctx
    print(f"[agent-endpoint] Set current_agent_ctx (idx={connection_idx})")

    async with agent_start_lock:
        if a_proc.conn_a is None:
            await start_agent_proc()

    conn_a = a_proc.conn_a
    if conn_a is None:
        # On reconnect, initial conn should not wipe global ctx on cleanup
        if entrypoint_module.current_agent_ctx is ctx:
            entrypoint_module.current_agent_ctx = None
        await ctx.send_message(
            orjson.dumps({
                "type": "agent_error",
                "error": "Agent process failed to start",
                "fatal": True
            }).decode()
        )
        return "Agent not available"

    connection_idx += 1

    try:
        while True:
            msg = await receive_json(ctx.receive)
            msg_type = msg.get("type", "unknown")
            tx_id = msg.get("tx_id", "none")
            print(f"[agent-endpoint] Received from browser: {msg_type} tx_id={tx_id}")
            await conn_a.send(msg)
            print(f"[agent-endpoint] Forwarded to agent: {msg_type} tx_id={tx_id}")
    finally:
        print(f"[agent-endpoint] Connection closed, clearing current_agent_ctx")
        if entrypoint_module.current_agent_ctx is ctx:
            entrypoint_module.current_agent_ctx = None

    return "Ok"
