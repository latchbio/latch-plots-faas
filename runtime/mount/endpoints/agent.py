import asyncio

import orjson
from latch_asgi.context.websocket import Context, HandlerResult
from latch_asgi.framework.websocket import (
    receive_json,
)
from latch_o11y.o11y import trace_app_function_with_span
from opentelemetry.trace import Span

import runtime.mount.entrypoint as entrypoint_module

from ..entrypoint import (
    a_proc,
    pod_id,
    pod_session_id,
    start_agent_proc,
    start_headless_browser,
)

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

    entrypoint_module.current_agent_ctx = ctx

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

            msg_type = msg.get("type")
            if msg_type == "init":
                local_storage = msg.get("local_storage")
                if isinstance(local_storage, dict):
                    local_storage.setdefault("plots.is_agent_controlled", "yes")
                    entrypoint_module.latest_local_storage = local_storage
                    print("[agent_ws] Stored latest auth state for headless browser hydration")

                if entrypoint_module.latest_local_storage is None:
                    print("[agent_ws] No auth state found, skipping headless browser start")
                    continue

                notebook_id = msg.get("notebook_id")
                try:
                    await start_headless_browser(
                        notebook_id,
                        local_storage=entrypoint_module.latest_local_storage,
                    )
                except Exception as e:
                    print(f"[agent_ws] Failed to start headless browser: {e}")

            await conn_a.send(msg)
    finally:
        if entrypoint_module.current_agent_ctx is ctx:
            entrypoint_module.current_agent_ctx = None

    return "Ok"
