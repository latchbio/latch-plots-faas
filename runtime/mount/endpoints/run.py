import re
import secrets
import signal
from contextlib import suppress
from dataclasses import dataclass

import aiohttp
import jwt
import orjson
from latch_asgi.context.websocket import Context, HandlerResult
from latch_asgi.framework.websocket import (
    WebsocketBadMessage,
    WebsocketConnectionClosedError,
    receive_json,
)
from latch_data_validation.data_validation import validate
from latch_o11y.o11y import trace_app_function_with_span
from opentelemetry.trace import Span

from ..entrypoint import (
    cell_last_run_outputs,
    cell_sequencers,
    cell_status,
    k_proc,
    plots_ctx_manager,
    pod_id,
    pod_session_id,
    ready_ev,
)
from ..utils import gql_query


@dataclass(frozen=True)
class MAuth:
    token: str
    notebook_id: str


@dataclass(frozen=True)
class PlotsSignerHasNotebookAccessData:
    plotsSignerHasNotebookAccess: bool


@dataclass(frozen=True)
class PlotsSignerHasNotebookAccessResp:
    data: PlotsSignerHasNotebookAccessData


async def raise_for_status(x: aiohttp.ClientResponse) -> None:
    if 200 <= x.status <= 299:
        return

    body = await x.read()

    try:
        data = orjson.loads(body)
        print(data)
    except orjson.JSONDecodeError:
        with suppress(UnicodeDecodeError):
            text = body.decode()
            print(text)

    x.raise_for_status()


connection_idx = 0

auth_header_regex = re.compile(
    r"""
    ^(
        Bearer \s+ (?P<oauth_token>.*) |
        Latch-Session-Token \s+ (?P<session_token>.*) |
    )$
    """,
    re.IGNORECASE | re.VERBOSE,
)


@trace_app_function_with_span
async def run(s: Span, ctx: Context) -> HandlerResult:
    global connection_idx, session_owner

    sess_hash = secrets.token_hex(32)
    await ctx.accept_connection()

    auth_msg = await ctx.receive_message(MAuth)

    notebook_id = auth_msg.notebook_id

    s.set_attributes(
        {
            "notebook_id": notebook_id,
            "pod_id": pod_id,
            "pod_session_id": pod_session_id,
            "sess_hash": sess_hash,
            "connection_idx": connection_idx,
        }
    )

    auth_header_regex_match = auth_header_regex.match(auth_msg.token)
    auth0_sub: str | None = None
    picture_url: str | None = None
    name: str | None = None

    if auth_header_regex_match is not None:
        oauth_token = auth_header_regex_match.group("oauth_token")
        session_token = auth_header_regex_match.group("session_token")

        if oauth_token is not None:
            # todo(rteqs): expose functionality in latch_asgi
            auth_data = jwt.decode(oauth_token, options={"verify_signature": False})
            auth0_sub = auth_data.get("sub")
            picture_url = auth_data.get("picture")
            name = auth_data.get("name")

            if auth0_sub is not None and name is not None:
                s.set_attribute("auth0_sub", auth0_sub)
                s.set_attribute("name", name)

        elif session_token is not None:
            s.set_attribute("session_token", session_token)

    data_q = await gql_query(
        auth=auth_msg.token,
        query="""
            query plotsSignerHasNotebookAccess($notebookId: BigInt!) {
                plotsSignerHasNotebookAccess(argNotebookId: $notebookId)
            }
        """,
        variables={"notebookId": notebook_id},
    )

    if "errors" in data_q:
        raise WebsocketBadMessage("failed to authenticate user")
    data = validate(data_q, PlotsSignerHasNotebookAccessResp)

    if not data.data.plotsSignerHasNotebookAccess:
        raise WebsocketBadMessage(f"Signer cannot access notebook {notebook_id}")

    conn_k = k_proc.conn_k
    assert conn_k is not None
    await plots_ctx_manager.add_context(
        sess_hash, ctx, connection_idx, auth0_sub, picture_url, name
    )

    await ready_ev.wait()

    await ctx.send_message(
        orjson.dumps(
            {
                "type": "ready",
                "connection_idx": connection_idx,
                "cell_status": cell_status,
                "cell_sequencers": cell_sequencers,
                "cell_outputs": cell_last_run_outputs,
            }
        ).decode()
    )

    connection_idx += 1

    try:
        while True:
            msg = await receive_json(ctx.receive)
            session_owner = plots_ctx_manager.session_owner
            is_session_owner = session_owner in {auth0_sub, connection_idx}

            if msg["type"] == "dispose_cell":
                cell_id = msg["cell_id"]

                if cell_status.get(cell_id) == "running" and k_proc.proc is not None:
                    k_proc.proc.send_signal(signal=signal.SIGINT)

                cell_status.pop(cell_id, None)
                cell_last_run_outputs.pop(cell_id, None)
                cell_sequencers.pop(cell_id, None)

            if msg["type"] == "stop_cell" and k_proc.proc is not None:
                k_proc.proc.send_signal(signal=signal.SIGINT)
                continue

            if msg["type"] == "run_cell" and not is_session_owner:
                await ctx.send_message(
                    orjson.dumps(
                        {
                            "type": "error",
                            "data": {"message": "user is not session owner"},
                        }
                    ).decode()
                )
                continue

            await conn_k.send(msg)
    except WebsocketConnectionClosedError:
        ...
    finally:
        await plots_ctx_manager.delete_context(sess_hash)

    return "Ok"
