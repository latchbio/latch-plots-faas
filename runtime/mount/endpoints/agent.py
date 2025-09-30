import asyncio
import os
import re
from dataclasses import dataclass

import jwt
import orjson
from latch_asgi.context.websocket import Context, HandlerResult
from latch_asgi.framework.websocket import (
    WebsocketBadMessage,
    receive_json,
)
from latch_data_validation.data_validation import validate
from latch_o11y.o11y import trace_app_function_with_span
from opentelemetry.trace import Span

from ..entrypoint import a_proc, pod_id, pod_session_id
from ..utils import gql_query

DISABLE_AUTH = os.environ.get("DISABLE_AUTH") == "1"


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
async def agent(s: Span, ctx: Context) -> HandlerResult:
    global connection_idx

    await ctx.accept_connection()

    auth_msg = await ctx.receive_message(MAuth)

    notebook_id = auth_msg.notebook_id

    s.set_attributes({
        "notebook_id": notebook_id,
        "pod_id": pod_id,
        "pod_session_id": pod_session_id,
        "connection_idx": connection_idx,
    })

    auth_header_regex_match = auth_header_regex.match(auth_msg.token)
    auth0_sub: str | None = None

    if auth_header_regex_match is not None:
        oauth_token = auth_header_regex_match.group("oauth_token")
        session_token = auth_header_regex_match.group("session_token")

        if oauth_token is not None:
            auth_data = jwt.decode(oauth_token, options={"verify_signature": False})
            auth0_sub = auth_data.get("sub")

            if auth0_sub is not None:
                s.set_attribute("auth0_sub", auth0_sub)

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

    conn_a = a_proc.conn_a
    if conn_a is None:
        await ctx.send_message(
            orjson.dumps({
                "type": "agent_error",
                "error": "Agent process not available",
                "fatal": True
            }).decode()
        )
        return "Agent not available"

    await conn_a.send({
        "type": "init",
        "context": f"You are assisting with notebook {notebook_id}."
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
