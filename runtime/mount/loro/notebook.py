import asyncio
import contextlib
import os
import random
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from aiohttp import (
    ClientConnectionError,
    ClientResponseError,
    ClientSession,
    ServerDisconnectedError,
)
from latch_data_validation.data_validation import validate
from yarl import URL

from loro import LoroDoc

latch_p = Path(os.environ.get("LATCH_SANDBOX_ROOT", "/root/.latch"))
sdk_token_path = latch_p / "token"
if sdk_token_path.exists():
    sdk_token = sdk_token_path.read_text()
else:
    sdk_token = ""
auth_token_sdk = f"Latch-SDK-Token {sdk_token}"

sess: ClientSession | None = None


def get_global_http_sess() -> ClientSession:
    global sess
    if sess is None:
        sess = ClientSession()

    return sess


def _retry_delay(attempt: int) -> float:
    delay = min(0.5 * (2**attempt), 30.0)
    jitter = delay * 0.5 * random.random()
    return delay + jitter


retryable_status_codes = frozenset({429, 500, 502, 503, 504})


def _should_retry(exc: Exception) -> bool:
    if isinstance(exc, ClientResponseError):
        return exc.status in retryable_status_codes
    return isinstance(exc, (ServerDisconnectedError, ClientConnectionError, OSError))


async def gql_query(
    query: str, variables: dict[str, Any], auth: str, *, max_retries: int = 5
) -> Any:
    sess = get_global_http_sess()
    url = URL("https://vacuole.latch.bio") / "graphql"
    headers = {
        "Content-Type": "application/json",
        "User-Agent": "LatchPlotsFaas/0.2.0",
        "Authorization": auth,
    }

    for attempt in range(max_retries + 1):
        try:
            async with sess.post(
                url, headers=headers, json={"query": query, "variables": variables}
            ) as resp:
                resp.raise_for_status()

                res = await resp.json()
                if "errors" in res:
                    raise RuntimeError(f"graphql error: {res}")
                return res

        except (
            ServerDisconnectedError,
            ClientConnectionError,
            ClientResponseError,
            OSError,
        ) as e:
            if _should_retry(e) and attempt < max_retries:
                delay = _retry_delay(attempt)

                if isinstance(e, ClientResponseError) and e.headers is not None:
                    retry_after = e.headers.get("Retry-After")
                    if retry_after is not None:
                        with contextlib.suppress(ValueError):
                            delay = max(delay, float(retry_after))

                print(
                    f"[gql_query] {type(e).__name__} on attempt"
                    f" {attempt + 1}/{max_retries + 1}, retrying in {delay:.1f}s: {e}"
                )
                await asyncio.sleep(delay)
                continue
            raise

    raise AssertionError("unreachable")


@dataclass
class NotebookCrdtUpdates:
    data: bytes


@dataclass
class GetNotebookCrdtUpdateNodes:
    nodes: list[NotebookCrdtUpdates]


@dataclass
class GetNotebookCrdtUpdatesData:
    plotNotebookCrdtUpdates: GetNotebookCrdtUpdateNodes


@dataclass
class GetNotebookCrdtUpdatesRes:
    data: GetNotebookCrdtUpdatesData


async def get_notebook_crdt_updates(
    notebook_id: str, latest_update_id: str | None
) -> list[bytes]:
    try:
        gql_res = await gql_query(
            query="""
                query GetPlotNotebookCrdtUpdates(
                    $notebookId: BigInt!
                    $latestUpdateId: BigInt
                ) {
                    plotNotebookCrdtUpdates(
                        filter: {
                            notebookId: { equalTo: $notebookId }
                            id: { greaterThan: $latestUpdateId }
                            loroVersion: { equalTo: "1.0" }
                        }
                    ) {
                        nodes {
                            data
                        }
                    }
                }
            """,
            variables={"notebookId": notebook_id, "latestUpdateId": latest_update_id},
            auth=auth_token_sdk,
        )
        updates = validate(gql_res, GetNotebookCrdtUpdatesRes)

        return [upd.data for upd in updates.plotNotebookCrdtUpdates.nodes]

    except Exception:
        # todo(rteqs): proper error handling
        traceback.print_exc()

    return []


@dataclass
class NotebookCheckpointInfo:
    data: bytes
    latestUpdateId: str


@dataclass
class GetPlotNotebookCheckpointNodes:
    nodes: list[NotebookCheckpointInfo]


@dataclass
class PlotNotebookInfo:
    data: bytes | None


@dataclass
class GetPlotNotebookCheckpointData:
    plotNotebookInfo: PlotNotebookInfo
    plotNotebookCheckpointInfos: GetPlotNotebookCheckpointNodes


@dataclass
class GetPlotNotebookCheckpointRes:
    data: GetPlotNotebookCheckpointData


# todo(rteqs): return type cringe
async def get_latest_checkpoint(notebook_id: str) -> tuple[bytes | None, str | None]:
    try:
        gql_res = await gql_query(
            query="""
                query GetPlotNotebookCheckpoint($notebookId: BigInt!) {
                    plotNotebookInfo(id: $notebookId) {
                        data
                    }
                    plotNotebookCheckpointInfos(
                        filter: {
                            notebookId: { equalTo: $notebookId }
                            loroVersion: { equalTo: "1.0" }
                        }
                        orderBy: [LATEST_UPDATE_ID_DESC, CREATION_TIME_DESC]
                        first: 1
                    ) {
                        nodes {
                            data
                            latestUpdateId
                        }
                    }
                }
            """,
            variables={"notebookId": notebook_id},
            auth=auth_token_sdk,
        )
        cp = validate(gql_res, GetPlotNotebookCheckpointRes)
        nodes = cp.plotNotebookCheckpointInfos.nodes

        if len(nodes) == 0:
            return (cp.plotNotebookInfo.data, None)

        n = cp.plotNotebookCheckpointInfos.nodes[0]
        return (n.data, n.latestUpdateId)

    except Exception:
        # todo(rteqs): proper error handling
        traceback.print_exc()
        raise


async def main() -> None:
    notebook_id = "52793"
    cp = await get_latest_checkpoint(notebook_id=notebook_id)

    crdt_updates = await get_notebook_crdt_updates(
        notebook_id=notebook_id, latest_update_id=cp[1]
    )

    batch = ([] if cp[0] is None else [cp[0]]) + crdt_updates  # todo(rteqs): cleanup

    doc = LoroDoc()
    doc.import_batch(batch)

    print(doc)
