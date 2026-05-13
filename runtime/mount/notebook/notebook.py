import base64
import contextlib
import json
import re
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from latch_data_validation.data_validation import validate
from loro import (
    ContainerID,
    ContainerType,
    ExportMode,
    LoroDoc,
    LoroMap,
    LoroMovableList,
    LoroText,
)
from utils import auth_token_sdk, gql_query


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
    notebook_id: str, latest_update_id: str | None
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
class PlotNotebookMetadataInfo:
    metadata: str | None


@dataclass(frozen=True)
class GetPlotNotebookMetadataData:
    plotNotebookInfo: PlotNotebookMetadataInfo | None


@dataclass(frozen=True)
class GetPlotNotebookMetadataRes:
    data: GetPlotNotebookMetadataData


@dataclass(frozen=True)
class CreatedPlotTransformInfo:
    id: str


@dataclass(frozen=True)
class CreatePlotTransformInfoPayload:
    plotTransformInfo: CreatedPlotTransformInfo


@dataclass(frozen=True)
class CreatePlotTransformInfoData:
    createPlotTransformInfo: CreatePlotTransformInfoPayload


@dataclass(frozen=True)
class CreatePlotTransformInfoRes:
    data: CreatePlotTransformInfoData


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


async def get_notebook_doc(notebook_id: str) -> NotebookCrdt:
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


class Notebook:
    def __init__(
        self,
        notebook_id: str,
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
    async def create(cls, notebook_id: str) -> "Notebook":
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
    ) -> tuple[str, str]:
        await self.fetch_updates()
        cell: LoroMap = self.cells.insert_container(pos, LoroMap())  # type: ignore

        cell.insert("cellType", "code")
        cell.insert("language", "python")

        source: LoroText = cell.insert_container("source", LoroText())  # type: ignore
        if code is not None:
            source.insert(0, code)

        cell_id = loro_ts_container_id(cell.id)
        try:
            gql_res = await gql_query(
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
                    "loroCellId": cell_id,
                    "parentNotebookId": self.notebook_id,
                },
                auth=auth_token_sdk,
            )
            res = validate(gql_res, CreatePlotTransformInfoRes)
            tf_id = res.data.createPlotTransformInfo.plotTransformInfo.id

        except Exception:
            # todo(rteqs): proper error handling
            traceback.print_exc()
            raise

        await self.export_updates()
        return cell_id, tf_id

    async def create_markdown_cell(self, pos: int, content: str | None) -> str:
        await self.fetch_updates()
        cell: LoroMap = self.cells.insert_container(pos, LoroMap())  # type: ignore

        cell.insert("cellType", "markdown")

        source: LoroText = cell.insert_container("source", LoroText())  # type: ignore
        if content is not None:
            source.insert(0, content)

        await self.export_updates()

        return loro_ts_container_id(cell.id)

    async def create_tab_marker_cell(self, pos: int, name: str) -> str:
        await self.fetch_updates()
        cell: LoroMap = self.cells.insert_container(pos, LoroMap())  # type: ignore

        cell.insert("cellType", "tabMarker")

        source: LoroText = cell.insert_container("displayName", LoroText())  # type: ignore
        if name is not None:
            source.insert(0, name)

        await self.export_updates()

        return loro_ts_container_id(cell.id)

    async def rename_tab(self, cell_id: str, new_name: str) -> None:
        await self.fetch_updates()

        if cell_id == "DEFAULT":
            gql_res = await gql_query(
                query="""
                    query GetPlotNotebookMetadata($notebookId: BigInt!) {
                        plotNotebookInfo(id: $notebookId) {
                            metadata
                        }
                    }
                """,
                variables={"notebookId": self.notebook_id},
                auth=auth_token_sdk,
            )
            res = validate(gql_res, GetPlotNotebookMetadataRes)

            info = res.data.plotNotebookInfo
            raw = info.metadata if info is not None else None
            meta: dict[str, Any] = json.loads(raw) if raw else {}
            meta["defaultTabName"] = new_name

            await gql_query(
                query="""
                    mutation PlotNotebookUpdateMetadata(
                        $id: BigInt!
                        $metadata: String!
                    ) {
                        updatePlotNotebookInfo(
                            input: { id: $id, patch: { metadata: $metadata } }
                        ) {
                            clientMutationId
                        }
                    }
                """,
                variables={"id": self.notebook_id, "metadata": json.dumps(meta)},
                auth=auth_token_sdk,
            )
            return

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
        await self.fetch_updates()
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
        await self.fetch_updates()
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
        await self.fetch_updates()
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

    async def rename_notebook(self, name: str) -> None:
        await gql_query(
            query="""
                mutation PlotsRenameNotebook($id: BigInt!, $displayName: String!) {
                    updatePlotNotebookInfo(
                        input: { id: $id, patch: { displayName: $displayName } }
                    ) {
                        clientMutationId
                    }
                }
            """,
            variables={"id": self.notebook_id, "displayName": name},
            auth=auth_token_sdk,
        )

    async def restore_checkpoint(self, template_version_id: str) -> None:
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
        await self.fetch_updates()


# todo(rteqs): remove duplication with entrypoint.py
latch_p = Path("/root/.latch")
notebook_id_file = latch_p / "notebook-id"
notebook_id: str | None = None
with contextlib.suppress(Exception):
    notebook_id = (
        notebook_id_file.read_text().strip() if notebook_id_file.exists() else None
    )

notebook: Notebook | None = None


async def get_notebook() -> Notebook:
    global notebook

    assert notebook_id is not None

    if notebook is None:
        notebook = await Notebook.create(notebook_id=notebook_id)

    return notebook
