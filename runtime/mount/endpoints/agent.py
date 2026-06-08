import asyncio

import orjson
from latch_asgi.context.websocket import Context, HandlerResult
from latch_asgi.framework.websocket import receive_json
from latch_o11y.o11y import trace_app_function_with_span
from opentelemetry.trace import Span

import runtime.mount.entrypoint as entrypoint_module

from ..entrypoint import (
    a_proc,
    handle_user_disconnect_fallback,
    mark_action_handled,
    pod_id,
    pod_session_id,
    start_agent_proc,
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

    async with agent_start_lock:
        if a_proc.msg_io is None:
            await start_agent_proc()

    conn_a = a_proc.msg_io
    if conn_a is None:
        await ctx.send_message(
            orjson.dumps({
                "type": "agent_error",
                "error": "Agent process failed to start",
                "fatal": True,
            }).decode()
        )
        return "Agent not available"

    connection_idx += 1

    try:
        while True:
            msg = await receive_json(ctx.receive)

            msg_type = msg.get("type")
            if msg_type == "init":
                entrypoint_module.user_agent_ctx = ctx

                notebook_id = msg.get("notebook_id")
                if notebook_id is None:
                    raise ValueError("notebook_id required in init message")

            if msg_type == "agent_action_response":
                tx_id = msg.get("tx_id")
                if tx_id is not None:
                    mark_action_handled(tx_id)

            await conn_a.send(msg)
    finally:
        entrypoint_module.user_agent_ctx = None
        # note(aidan): not sure if this is preferable to the agent figuring it out
        await handle_user_disconnect_fallback()

    return "Ok"
