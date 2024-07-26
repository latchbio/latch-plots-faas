import json
import secrets
from contextlib import suppress
from dataclasses import dataclass

import aiohttp
import orjson
from latch_asgi.context.websocket import Context, HandlerResult
from latch_asgi.framework.websocket import (
    WebsocketBadMessage,
    WebsocketConnectionClosedError,
    receive_json,
)
from latch_data_validation.data_validation import validate

from ..entrypoint import (
    cell_last_run_outputs,
    cell_sequencers,
    cell_status,
    contexts,
    gql_query,
    k_proc,
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
active_cell: str | None = None


async def run(ctx: Context) -> HandlerResult:
    global connection_idx, active_cell

    sess_hash = secrets.token_hex(32)
    await ctx.accept_connection()

    auth_msg = await ctx.receive_message(MAuth)

    notebook_id = auth_msg.notebook_id

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
    connection_idx += 1

    if active_cell is not None and cell_status[active_cell] == "running":
        await ctx.send_message({"type": "cell_running", "cell_id": active_cell})

    try:
        while True:
            msg = await receive_json(ctx.receive)

            if msg["type"] == "dispose_cell":
                cell_id = msg["cell_id"]
                cell_status.pop(cell_id, None)
                cell_last_run_outputs.pop(cell_id, None)
                cell_sequencers.pop(cell_id, None)

            if msg["type"] == "run_cell":
                cell_status[msg["cell_id"]] = "running"
                active_cell = msg["cell_id"]

            await conn_k.send(msg)
    except WebsocketConnectionClosedError:
        ...
    finally:
        del contexts[sess_hash]

    return "Ok"
