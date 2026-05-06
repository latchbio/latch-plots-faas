import traceback
from dataclasses import dataclass

from latch_data_validation.data_validation import validate
from utils import auth_token_sdk

from loro import LoroDoc
from runtime.mount.utils import gql_query


@dataclass
class NotebookCrdtUpdates:
    data: bytes


@dataclass
class GetNotebookCrdtUpdateNodes:
    nodes: list[NotebookCrdtUpdates]


@dataclass
class GetNotebookCrdtUpdatesRes:
    plotNotebookCrdtUpdates: GetNotebookCrdtUpdateNodes


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
class GetPlotNotebookCheckpointRes:
    plotNotebookInfo: PlotNotebookInfo
    plotNotebookCheckpointInfos: GetPlotNotebookCheckpointNodes


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
