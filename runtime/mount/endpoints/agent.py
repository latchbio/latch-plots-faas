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
    pod_id,
    pod_session_id,
    start_agent_proc,
    start_headless_browser,
)

connection_idx = 0
agent_start_lock = asyncio.Lock()

# Connection role is determined by the init message:
# - "user" browser sends init WITH local_storage (triggers headless browser start)
# - "action_handler" (headless browser) sends init WITHOUT local_storage
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

    # Don't set context yet - wait for init message to determine role
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
                    # User browser: sends init WITH local_storage (includes auth data)
                    # This triggers headless browser startup
                    connection_role = "user"
                    entrypoint_module.user_agent_ctx = ctx
                    entrypoint_module.latest_local_storage = local_storage

                    notebook_id = msg.get("notebook_id")
                    if notebook_id is not None and entrypoint_module.latest_local_storage is not None:
                        try:
                            await start_headless_browser(
                                notebook_id,
                                local_storage=entrypoint_module.latest_local_storage,
                            )
                        except Exception as e:
                            print(f"[agent_ws] Failed to start headless browser: {e}")
                else:
                    # Headless browser: sends init WITHOUT local_storage
                    # This is the action handler that executes agent actions
                    connection_role = "action_handler"
                    entrypoint_module.action_handler_ctx = ctx
                    action_handler_ready_ev.set()

            await conn_a.send(msg)
    finally:
        # Only clear the context that this connection owns
        if connection_role == "user" and entrypoint_module.user_agent_ctx is ctx:
            entrypoint_module.user_agent_ctx = None
        elif connection_role == "action_handler" and entrypoint_module.action_handler_ctx is ctx:
            entrypoint_module.action_handler_ctx = None
            action_handler_ready_ev.clear()

    return "Ok"
