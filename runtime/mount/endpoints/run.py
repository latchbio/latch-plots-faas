import json
import re
import secrets
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
    contexts,
    gql_query,
    k_proc,
    pod_id,
    pod_session_id,
    ready_ev,
)


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
    global connection_idx

    connection_idx += 1

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
    if auth_header_regex_match is not None:
        oauth_token = auth_header_regex_match.group("oauth_token")
        session_token = auth_header_regex_match.group("session_token")

        if oauth_token is not None:
            auth_data = jwt.decode(oauth_token, options={"verify_signature": False})
            s.set_attribute("auth0_sub", auth_data.get("sub"))
            s.set_attribute("name", auth_data.get("name"))
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
    contexts[sess_hash] = ctx

    await ready_ev.wait()
    await ctx.send_message(
        json.dumps(
            {
                "type": "ready",
                "connection_idx": connection_idx,
                "cell_status": cell_status,
                "cell_sequencers": cell_sequencers,
                "cell_outputs": cell_last_run_outputs,
            },
            allow_nan=False,
        )
    )

    try:
        while True:
            msg = await receive_json(ctx.receive)

            if msg["type"] == "dispose_cell":
                cell_id = msg["cell_id"]
                cell_status.pop(cell_id, None)
                cell_last_run_outputs.pop(cell_id, None)
                cell_sequencers.pop(cell_id, None)

            await conn_k.send(msg)
    except WebsocketConnectionClosedError:
        ...
    finally:
        del contexts[sess_hash]

    return "Ok"
