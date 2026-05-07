import asyncio
import base64
import contextlib
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

from loro import ContainerID, ExportMode, LoroDoc, LoroMap, LoroMovableList, LoroText

latch_p = Path("../../../scratch-local")
sdk_token_path = latch_p / "token"
if sdk_token_path.exists():
    sdk_token = sdk_token_path.read_text().strip()
else:
    raise Exception("SDK Token not found")
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


@dataclass(frozen=True)
class NotebookCrdtUpdates:
    id: str
    data: str


@dataclass(frozen=True)
class GetNotebookCrdtUpdateNodes:
    nodes: list[NotebookCrdtUpdates]


@dataclass(frozen=True)
class GetNotebookCrdtUpdatesData:
    plotNotebookCrdtUpdates: GetNotebookCrdtUpdateNodes


@dataclass(frozen=True)
class GetNotebookCrdtUpdatesRes:
    data: GetNotebookCrdtUpdatesData


@dataclass(frozen=True, kw_only=True)
class CrdtUpdates:
    updates: list[bytes]
    latest_update_id: int | None


async def get_notebook_crdt_updates(
    notebook_id: int, latest_update_id: int | None
) -> CrdtUpdates:
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
                            id
                            data
                        }
                    }
                }
            """,
            variables={"notebookId": notebook_id, "latestUpdateId": latest_update_id},
            auth=auth_token_sdk,
        )
        res = validate(gql_res, GetNotebookCrdtUpdatesRes)
        nodes = res.data.plotNotebookCrdtUpdates.nodes

        return CrdtUpdates(
            updates=[base64.b64decode(upd.data) for upd in nodes],
            latest_update_id=int(nodes[-1].id) if len(nodes) > 0 else latest_update_id,
        )

    except Exception:
        # todo(rteqs): proper error handling
        traceback.print_exc()
        raise


@dataclass(frozen=True)
class NotebookCheckpointInfo:
    data: str
    latestUpdateId: str


@dataclass(frozen=True)
class GetPlotNotebookCheckpointNodes:
    nodes: list[NotebookCheckpointInfo]


@dataclass(frozen=True)
class PlotNotebookInfo:
    data: str | None


@dataclass(frozen=True)
class GetPlotNotebookCheckpointData:
    plotNotebookInfo: PlotNotebookInfo
    plotNotebookCheckpointInfos: GetPlotNotebookCheckpointNodes


@dataclass(frozen=True)
class GetPlotNotebookCheckpointRes:
    data: GetPlotNotebookCheckpointData


@dataclass(frozen=True, kw_only=True)
class NotebookCrdt:
    loro_doc: LoroDoc
    latest_update_id: int | None


async def get_notebook_doc(notebook_id: int) -> NotebookCrdt:
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
        res = validate(gql_res, GetPlotNotebookCheckpointRes)

        import_bytes: list[bytes] = []

        cp_nodes = res.data.plotNotebookCheckpointInfos.nodes
        latest_update_id: int | None = None
        if len(cp_nodes) == 0:
            snapshot = res.data.plotNotebookInfo.data
            if snapshot is not None:
                import_bytes.append(base64.b64decode(snapshot))
        else:
            cp = cp_nodes[0]
            import_bytes.append(base64.b64decode(cp.data))
            latest_update_id = int(cp.latestUpdateId)

        crdt_updates = await get_notebook_crdt_updates(
            notebook_id=notebook_id, latest_update_id=latest_update_id
        )
        import_bytes.extend(crdt_updates.updates)

        doc = LoroDoc()
        doc.import_batch(bytes=import_bytes)

        return NotebookCrdt(
            loro_doc=doc, latest_update_id=crdt_updates.latest_update_id
        )

    except Exception:
        # todo(rteqs): proper error handling
        traceback.print_exc()
        raise


class Notebook:
    def __init__(
        self, notebook_id: int, loro_doc: LoroDoc, latest_update_id: int | None
    ) -> None:
        self.notebook_id = notebook_id
        self.loro_doc = loro_doc
        self.latest_update_id = latest_update_id
        self.last_persisted_version = self.loro_doc.oplog_vv

    @classmethod
    async def create(cls, notebook_id: int) -> "Notebook":
        res = await get_notebook_doc(notebook_id)
        return cls(
            notebook_id=notebook_id,
            loro_doc=res.loro_doc,
            latest_update_id=res.latest_update_id,
        )

    async def fetch_updates(self) -> None:
        res = await get_notebook_crdt_updates(
            notebook_id=self.notebook_id, latest_update_id=self.latest_update_id
        )

        self.latest_update_id = res.latest_update_id
        self.loro_doc.import_batch(res.updates)

    async def export_updates(self) -> None:
        upd = self.loro_doc.export(mode=ExportMode.Updates(self.last_persisted_version))
        self.last_persisted_version = self.loro_doc.oplog_vv
        try:
            await gql_query(
                query="""
                    mutation CreatePlotNotebookCrdtUpdate(
                        $notebookId: BigInt!
                        $data: Base64EncodedBinary!
                    ) {
                        createPlotNotebookCrdtUpdate(
                            input: {
                                plotNotebookCrdtUpdate: {
                                    notebookId: $notebookId
                                    data: $data
                                    loroVersion: "1.0"
                                }
                            }
                        ) {
                            clientMutationId
                        }
                    }
                """,
                variables={
                    "notebookId": self.notebook_id,
                    "data": base64.b64encode(upd).decode(),
                },
                auth=auth_token_sdk,
            )

        except Exception:
            # todo(rteqs): proper error handling
            traceback.print_exc()
            raise

    @property
    def cells(self) -> LoroMovableList:
        return self.loro_doc.get_movable_list("cells")


def loro_ts_container_id(c: ContainerID) -> str:
    if isinstance(c, str):
        return c

    if isinstance(c, c.Root):
        return f"cid:root-{c.name}:{c.container_type}"

    if isinstance(c, c.Normal):
        return f"cid:{c.counter}@{c.peer}:{c.container_type}"

    # todo(rteqs): proper exception
    raise


async def create_code_cell(
    notebook: Notebook, pos: int, code: str | None, display_name: str
) -> None:
    cell: LoroMap = notebook.cells.insert_container(pos, LoroMap())  # type: ignore

    cell.insert("cellType", "code")
    cell.insert("language", "python")

    source = LoroText()
    if code is not None:
        source.insert(0, code)
    cell.insert_container("source", source)

    try:
        await gql_query(
            query="""
                mutation PlotsCreateTransform(
                    $ownerId: BigInt!
                    $displayName: String!
                    $parentNotebookId: BigInt!
                    $loroCellId: String!
                ) {
                    createPlotTransformInfo(
                        input: {
                            plotTransformInfo: {
                                ownerId: $ownerId
                                displayName: $displayName
                                parentNotebookId: $parentNotebookId
                                loroCellId: $loroCellId
                            }
                        }
                    ) {
                        clientMutationId
                        plotTransformInfo {
                            id
                        }
                    }
                }
            """,
            variables={
                # todo(rteqs): replace with account id from accountInfoCurrent query
                "ownerId": 22657,
                "displayName": display_name,
                "loroCellId": loro_ts_container_id(cell.id),
                "parentNotebookId": notebook.notebook_id,
            },
            auth=auth_token_sdk,
        )

    except Exception:
        # todo(rteqs): proper error handling
        traceback.print_exc()
        raise

    await notebook.export_updates()


# todo(rteqs): implement
async def get_account_info(): ...


async def main() -> None:
    notebook = await Notebook.create(52793)
    pos = len(notebook.cells)
    await create_code_cell(
        notebook=notebook, pos=pos, code="print(100)", display_name="Untitled Transform"
    )

    if sess is not None:
        await sess.close()


asyncio.run(main())
