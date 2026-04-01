import asyncio
from typing import Literal

import orjson
from latch_asgi.context.websocket import Context, HandlerResult
from latch_asgi.framework.websocket import receive_json
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
from ..utils import auth_token_sdk, gql_query, sdk_token

connection_idx = 0
agent_start_lock = asyncio.Lock()

ConnectionRole = Literal["user", "action_handler", "unknown"]


async def _resolve_headless_browser_bootstrap() -> tuple[str, dict[str, str]] | None:
    if pod_id is None:
        return None

    resp = await gql_query(
        auth=auth_token_sdk,
        query="""
            query AgentBootstrapPodInfo($podId: BigInt!) {
                podInfos(filter: { id: { equalTo: $podId } }, first: 1) {
                    nodes {
                        plotNotebookId
                        accountId
                    }
                }
            }
        """,
        variables={"podId": pod_id},
    )

    data = resp.get("data")
    if not isinstance(data, dict):
        return None

    pod_infos = data.get("podInfos")
    if not isinstance(pod_infos, dict):
        return None

    nodes = pod_infos.get("nodes")
    if not isinstance(nodes, list) or len(nodes) == 0:
        return None

    node = nodes[0]
    if not isinstance(node, dict):
        return None

    notebook_id = node.get("plotNotebookId")
    workspace_id = node.get("accountId")
    if notebook_id is None or workspace_id is None:
        return None

    local_storage = {
        "plots.is_agent_controlled": "yes",
        "viewAccountId": str(workspace_id),
        "latch.authData": orjson.dumps({
            "status": "done",
            "auth0Data": {
                "idToken": sdk_token.strip(),
                "idTokenPayload": {
                    "sub": "agent-session",
                    "latch.bio/tos_ok": "true",
                },
            },
        }).decode(),
    }
    return str(notebook_id), local_storage


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
        bootstrap = await _resolve_headless_browser_bootstrap()
        if bootstrap is not None:
            print("[agent endpoint] Starting headless browser")
            notebook_id, local_storage = bootstrap
            await start_headless_browser(notebook_id, local_storage=local_storage)
        else:
            raise RuntimeError("Failed to resolve headless browser bootstrap")
    except Exception as e:  # noqa: BLE001
        print(f"[agent endpoint] Failed to bootstrap headless browser: {e!s}")

    try:
        while True:
            msg = await receive_json(ctx.receive)

            msg_type = msg.get("type")
            if msg_type == "init":
                local_storage = msg.get("local_storage")

                if isinstance(local_storage, dict) and len(local_storage) > 0:
                    connection_role = "user"
                    entrypoint_module.user_agent_ctx = ctx

                    notebook_id = msg.get("notebook_id")
                    if notebook_id is None:
                        raise ValueError("notebook_id required in init message")

                    # This will be a noop if the browser is already running with the same local storage or if its an agent calling init
                    # If this is called by a user session and the local storage has changed then the browser will be restarted
                    await start_headless_browser(
                        notebook_id,
                        local_storage=local_storage,
                    )
                else:
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
            # note(aidan): not sure if this is preferable to the agent figuring it out
            await handle_user_disconnect_fallback()
        elif connection_role == "action_handler":
            entrypoint_module.action_handler_ctx = None
            action_handler_ready_ev.clear()
            await restart_headless_browser()

    return "Ok"
