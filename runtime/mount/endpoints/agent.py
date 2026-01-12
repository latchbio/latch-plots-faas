import asyncio
from typing import Literal

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
    action_handler_ready_ev,
    handle_user_disconnect_fallback,
    mark_action_handled,
    pod_id,
    pod_session_id,
    restart_headless_browser,
    start_agent_proc,
    start_headless_browser,
)

connection_idx = 0
agent_start_lock = asyncio.Lock()

ConnectionRole = Literal["user", "action_handler", "unknown"]


@trace_app_function_with_span
async def agent(s: Span, ctx: Context) -> HandlerResult:
    global connection_idx

    await ctx.accept_connection()

    s.set_attributes({
        "pod_id": pod_id,
        "pod_session_id": pod_session_id,
        "connection_idx": connection_idx,
    })

    connection_role: ConnectionRole = "unknown"

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

    connection_idx += 1

    try:
        while True:
            msg = await receive_json(ctx.receive)

            msg_type = msg.get("type")
            if msg_type == "init":
                local_storage = msg.get("local_storage")

                if isinstance(local_storage, dict) and len(local_storage) > 0:
                    connection_role = "user"
                    entrypoint_module.user_agent_ctx = ctx
                    entrypoint_module.latest_local_storage = local_storage

                    notebook_id = msg.get("notebook_id")
                    if notebook_id is None:
                        raise ValueError("notebook_id required in init message")

                    await start_headless_browser(
                        notebook_id,
                        local_storage=entrypoint_module.latest_local_storage,
                    )
                else:
                    # backend browser opening connection
                    connection_role = "action_handler"
                    entrypoint_module.action_handler_ctx = ctx
                    action_handler_ready_ev.set()

            if msg_type == "agent_action_response":
                tx_id = msg.get("tx_id")
                if tx_id is not None:
                    mark_action_handled(tx_id)

            await conn_a.send(msg)
    finally:
        if connection_role == "user":
            entrypoint_module.user_agent_ctx = None
            # Fallback any pending user browser actions to backend browser
            await handle_user_disconnect_fallback()
        elif connection_role == "action_handler":
            entrypoint_module.action_handler_ctx = None
            action_handler_ready_ev.clear()
            await restart_headless_browser()

    return "Ok"
