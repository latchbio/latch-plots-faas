import asyncio
import base64
import contextlib
import random
import re
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

from loro import (
    ContainerID,
    ContainerType,
    ExportMode,
    LoroDoc,
    LoroMap,
    LoroMovableList,
    LoroText,
)

# todo(rteqs): a bunch of the loro stuff we did here could be upstreamed

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
    latest_update_id: str | None


async def get_notebook_crdt_updates(
    notebook_id: int, latest_update_id: str | None
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
            latest_update_id=nodes[-1].id if len(nodes) > 0 else latest_update_id,
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
    ownerId: str
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
    latest_update_id: str | None
    owner_id: str


async def get_notebook_doc(notebook_id: int) -> NotebookCrdt:
    try:
        gql_res = await gql_query(
            query="""
                query GetPlotNotebookCheckpoint($notebookId: BigInt!) {
                    plotNotebookInfo(id: $notebookId) {
                        ownerId
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

        cp_nodes = res.data.plotNotebookCheckpointInfos.nodes
        latest_update_id: str | None = None
        doc = LoroDoc()

        if len(cp_nodes) == 0:
            snapshot = res.data.plotNotebookInfo.data
            if snapshot is not None:
                doc.import_(base64.b64decode(snapshot))
        else:
            cp = cp_nodes[0]
            doc.import_(base64.b64decode(cp.data))
            latest_update_id = cp.latestUpdateId

        crdt_updates = await get_notebook_crdt_updates(
            notebook_id=notebook_id, latest_update_id=latest_update_id
        )
        doc.import_batch(crdt_updates.updates)

        return NotebookCrdt(
            loro_doc=doc,
            latest_update_id=crdt_updates.latest_update_id,
            owner_id=res.data.plotNotebookInfo.ownerId,
        )

    except Exception:
        # todo(rteqs): proper error handling
        traceback.print_exc()
        raise


def loro_ts_container_id(c: ContainerID) -> str:
    if isinstance(c, str):
        return c

    if isinstance(c, c.Root):
        return f"cid:root-{c.name}:{c.container_type}"

    if isinstance(c, c.Normal):
        return f"cid:{c.counter}@{c.peer}:{c.container_type}"

    # todo(rteqs): proper exception
    raise


root_container_id_re = re.compile(
    r"""
       ^cid:               # prefix
       root-               # root marker
       (?P<name>.+)        # root name (any chars)
       :(?P<type>\w+)$     # container type
       """,
    re.VERBOSE,
)

normal_container_id_re = re.compile(
    r"""
       ^cid:               # prefix
       (?P<counter>\d+)    # counter
       @(?P<peer>\d+)      # peer id
       :(?P<type>\w+)$     # container type
       """,
    re.VERBOSE,
)

container_type_map = {
    "Text": ContainerType.Text(),
    "Map": ContainerType.Map(),
    "List": ContainerType.List(),
    "MovableList": ContainerType.MovableList(),
    "Tree": ContainerType.Tree(),
    "Counter": ContainerType.Counter(),
}


def parse_ts_container_id(s: str) -> ContainerID:
    root_match = root_container_id_re.match(s)
    if root_match is not None:
        return ContainerID.Root(
            name=root_match.group("name"),
            container_type=container_type_map[root_match.group("type")],
        )

    normal_match = normal_container_id_re.match(s)
    if normal_match is not None:
        return ContainerID.Normal(
            peer=int(normal_match.group("peer")),
            counter=int(normal_match.group("counter")),
            container_type=container_type_map[normal_match.group("type")],
        )

    raise
    # if s.contains


# class ContainerID:
#     class Root(ContainerID):
#         def __init__(self, name: str, container_type: ContainerType): ...
#         name: str
#         container_type: ContainerType

#     class Normal(ContainerID):
#         def __init__(self, peer: int, counter: int, container_type: ContainerType): ...
#         peer: int
#         counter: int
#         container_type: ContainerType


class Notebook:
    def __init__(
        self,
        notebook_id: int,
        loro_doc: LoroDoc,
        latest_update_id: str | None,
        owner_id: str,
    ) -> None:
        self.notebook_id = notebook_id
        self.loro_doc = loro_doc
        self.latest_update_id = latest_update_id
        self.last_persisted_version = self.loro_doc.oplog_vv
        self.owner_id = owner_id

    @classmethod
    async def create(cls, notebook_id: int) -> "Notebook":
        res = await get_notebook_doc(notebook_id)
        return cls(
            notebook_id=notebook_id,
            loro_doc=res.loro_doc,
            latest_update_id=res.latest_update_id,
            owner_id=res.owner_id,
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

    async def create_code_cell(
        self, pos: int, code: str | None, display_name: str
    ) -> None:
        cell: LoroMap = self.cells.insert_container(pos, LoroMap())  # type: ignore

        cell.insert("cellType", "code")
        cell.insert("language", "python")

        source: LoroText = cell.insert_container("source", LoroText())  # type: ignore
        if code is not None:
            source.insert(0, code)

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
                    "ownerId": self.owner_id,
                    "displayName": display_name,
                    "loroCellId": loro_ts_container_id(cell.id),
                    "parentNotebookId": self.notebook_id,
                },
                auth=auth_token_sdk,
            )

        except Exception:
            # todo(rteqs): proper error handling
            traceback.print_exc()
            raise

        await self.export_updates()

    async def create_markdown_cell(self, pos: int, content: str | None) -> None:
        cell: LoroMap = self.cells.insert_container(pos, LoroMap())  # type: ignore

        cell.insert("cellType", "markdown")

        source: LoroText = cell.insert_container("source", LoroText())  # type: ignore
        if content is not None:
            source.insert(0, content)

        await self.export_updates()

    async def create_tab_marker_cell(self, pos: int, name: str) -> None:
        cell: LoroMap = self.cells.insert_container(pos, LoroMap())  # type: ignore

        cell.insert("cellType", "tabMarker")

        source: LoroText = cell.insert_container("displayName", LoroText())  # type: ignore
        if name is not None:
            source.insert(0, name)

        await self.export_updates()

    async def rename_tab(self, cell_id: str, new_name: str) -> None:
        cell = self.loro_doc.get_container(id=parse_ts_container_id(cell_id))
        if cell is None:
            raise

        assert isinstance(cell, LoroMap)
        assert cell.get("cellType") == "tabMarker"

        source = cell.get("source")
        if source is None:
            raise

        assert isinstance(source, LoroText)
        source.update(new_name)

        await self.export_updates()

    async def edit_cell(self, cell_id: str, new_code: str) -> None:
        cell = self.loro_doc.get_container(id=parse_ts_container_id(cell_id))
        if cell is None:
            raise

        assert isinstance(cell, LoroMap)
        assert cell.get("cellType") == "code"

        source = cell.get("source")
        if source is None:
            raise

        assert isinstance(source, LoroText)
        source.update(new_code)

        await self.export_updates()

    async def delete_cell(self, cell_id: str) -> None:
        pos = -1
        for i, cid in enumerate(self.cells.to_vec()):
            if loro_ts_container_id(cid) == cell_id:
                pos = i
                break

        self.cells.delete(pos, 1)
        await self.export_updates()

        await gql_query(
            query="""
                mutation DeleteCell($loroCellId: String!) {
                    plotsTransformDeleteByLoroCellId(input: { loroCellId: $loroCellId }) {
                        clientMutationId
                    }
                    plotsDeleteByLoroCellId(input: { loroCellId: $loroCellId }) {
                        clientMutationId
                    }
                    plotCellValueViewerDeleteByLoroCellId(input: { loroCellId: $loroCellId }) {
                        clientMutationId
                    }
                }
            """,
            variables={"loroCellId": cell_id},
            auth=auth_token_sdk,
        )

    async def delete_all_cells(self) -> None:
        self.cells.clear()
        await self.export_updates()

        await gql_query(
            query="""
                mutation DeleteAllCells($notebookId: BigInt!) {
                    deleteAllCells(input: { argNotebookId: $notebookId }) {
                        clientMutationId
                    }
                }
            """,
            variables={"notebookId": self.notebook_id},
            auth=auth_token_sdk,
        )

    async def restore_checkpoint(self, template_version_id: str):
        await gql_query(
            query="""
                mutation PlotsRestoreNotebookToVersion($plotTemplateVersionId: BigInt!) {
                    plotsRestoreNotebookToVersion(
                        input: { argPlotTemplateVersionId: $plotTemplateVersionId }
                    ) {
                        clientMutationId
                    }
                }
            """,
            variables={"plotTemplateVersionId": template_version_id},
            auth=auth_token_sdk,
        )


async def main() -> None:
    notebook = await Notebook.create(52793)
    # print(notebook.owner_id)
    # print(
    #     notebook.loro_doc.get_container(
    #         id=parse_ts_container_id("cid:0@18108464718435593442:Map")
    #     ).get_deep_value()
    # )
    pos = 0
    for pos, cid in enumerate(notebook.cells.to_vec()):
        if loro_ts_container_id(cid) == "cid:0@18108464718435593442:Map":
            print("found")
            break
    print(pos)

    # cells = notebook.cells
    # for i in range(len(cells)):
    #     print(cells.get(i).id)
    # if c is None:
    #     continue
    # if c.is_container(c):
    #     print(c.container.id)

    # pos = len(notebook.cells)
    # print(notebook.loro_doc.get_deep_value())
    # await notebook.create_code_cell(pos=pos, code="print(100)", display_name="Test")
    # await notebook.create_markdown_cell(pos=pos, content="# Markdown Example")

    if sess is not None:
        await sess.close()


asyncio.run(main())
