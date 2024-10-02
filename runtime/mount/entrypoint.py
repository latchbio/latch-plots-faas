import asyncio
import contextlib
import os
import socket
import sys
import traceback
from asyncio.subprocess import Process
from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict, TypeVar

import orjson
from duckdb import DuckDBPyConnection
from latch_asgi.framework.websocket import WebsocketConnectionClosedError
from latch_data_validation.data_validation import validate

from .socketio import SocketIo
from .utils import PlotConfig, get_global_http_sess, gql_query

dir_p = Path(__file__).parent


from latch_asgi.context.websocket import Context

T = TypeVar("T")


sock, sock_k = socket.socketpair(family=socket.AF_UNIX)
sock.setblocking(False)
sock_k_fd = sock_k.detach()

ready_ev = asyncio.Event()

active_cell: str | None = None


class CellOutputs(TypedDict):
    outputs: list[str]
    dataframe_outputs: list[str]
    figure_outputs: list[str]


cell_status: dict[str, str] = {}
cell_sequencers: dict[str, int] = {}
cell_last_run_outputs: dict[str, CellOutputs] = {}

async_tasks: list[asyncio.Task] = []

contexts: dict[str, Context] = {}

latch_p = Path("/root/.latch")
sdk_token = (latch_p / "token").read_text()
auth_token_sdk = f"Latch-SDK-Token {sdk_token}"

pod_id = int((latch_p / "id").read_text())
pod_session_id = (latch_p / "session-id").read_text()


@dataclass
class KernelProc:
    conn_k: SocketIo | None = None
    proc: Process | None = None


k_proc = KernelProc()


# todo(maximsmol): typing
# todo(maximsmol): share types with kernel.py
class PaginationSettingsData(TypedDict):
    sort_settings: object | None
    row_filters: object | None
    selections: object | None


class ViewerCellData(TypedDict):
    source: list[str]  # todo(maximsmol): tuple[str, KeyType]
    pagination_settings: PaginationSettingsData


@dataclass(frozen=True)
class KernelState:
    widget_states: dict[str, str]
    cell_output_selections: dict[str, str]
    plot_data_selections: dict[str, str]
    viewer_cell_data: dict[str, ViewerCellData]
    plot_configs: dict[str, PlotConfig]


@dataclass(frozen=True)
class PlotsNotebookKernelState:
    plotsNotebookKernelState: KernelState


@dataclass(frozen=True)
class PlotsNotebookKernelStateResp:
    data: PlotsNotebookKernelState


conn: DuckDBPyConnection | None = None


async def try_send_message(ctx: Context, msg: str) -> None:
    with contextlib.suppress(WebsocketConnectionClosedError):
        await ctx.send_message(msg)


async def broadcast_message(msg: str) -> None:
    tasks = [
        asyncio.create_task(try_send_message(ctx, msg)) for ctx in contexts.values()
    ]
    await asyncio.gather(*tasks)


async def add_pod_event(*, auth: str, event_type: str) -> None:
    try:
        await gql_query(
            auth=auth,
            query="""
                mutation AddSessionEvent($eventType: String!, $podSessionId: BigInt!) {
                    createPodSessionEvent(
                        input: {podSessionEvent: {podSessionId: $podSessionId, eventType: $eventType}}
                    ) {
                        clientMutationId
                    }
                }

            """,
            variables={"eventType": event_type, "podSessionId": pod_session_id},
        )
    except Exception:
        traceback.print_exc()


async def handle_kernel_messages(conn_k: SocketIo, auth: str) -> None:
    global active_cell

    print("Starting kernel message listener")
    while True:
        msg = await conn_k.recv()
        print(">", msg)

        if msg["type"] == "ready":
            ready_ev.set()
            await add_pod_event(auth=auth, event_type="kernel_ready")
            continue

        if msg["type"] == "debug_state":
            msg["sess_hash"] = pod_session_id

        elif msg["type"] == "start_cell":
            cell_id = msg["cell_id"]
            active_cell = cell_id
            cell_sequencers[cell_id] = msg["run_sequencer"]
            cell_status[cell_id] = "running"

            await gql_query(
                auth=auth,
                query="""
                    mutation ClearCellMetadata($id: BigInt!) {
                        updatePlotTransformInfo(
                            input: { id: $id, patch: { logs: "", exception: null } }
                        ) {
                            clientMutationId
                        }
                    }
                """,
                variables={"id": msg["cell_id"]},
            )

            msg = {
                "type": msg["type"],
                "cell_id": msg["cell_id"],
                "run_sequencer": msg["run_sequencer"],
            }

        elif msg["type"] == "cell_result":
            cell_id = msg["cell_id"]
            cell_status[cell_id] = "ran" if "exception" not in msg else "error"
            cell_last_run_outputs[cell_id] = {
                "outputs": msg["outputs"],
                "dataframe_outputs": msg["dataframe_outputs"],
                "figure_outputs": msg["figure_outputs"],
            }

            exc = msg.get("exception")
            if exc is not None:
                exc = orjson.dumps({"string": exc[-9000:]}).decode()

            await gql_query(
                auth=auth,
                query="""
                    mutation UpdateCellResult($id: BigInt!, $exception: String) {
                        updatePlotTransformInfo(
                            input: { id: $id, patch: { exception: $exception } }
                        ) {
                            clientMutationId
                        }
                    }
                """,
                variables={"id": msg["cell_id"], "exception": exc},
            )

            msg = {
                "type": msg["type"],
                "cell_id": msg["cell_id"],
                "has_exception": exc is not None,
                "outputs": msg.get("outputs"),
                "dataframe_outputs": msg.get("dataframe_outputs"),
                "figure_outputs": msg.get("figure_outputs"),
            }

        elif msg["type"] == "cell_widgets":
            await gql_query(
                auth=auth,
                query="""
                    mutation UpdateCellWidgets($id: BigInt!, $state: String!) {
                        updatePlotTransformInfo(
                            input: { id: $id, patch: { widgetState: $state } }
                        ) {
                            clientMutationId
                        }
                    }
                """,
                variables={
                    "id": msg["cell_id"],
                    "state": orjson.dumps(msg["widget_state"]).decode(),
                },
            )

            msg = {
                "type": msg["type"],
                "cell_id": msg["cell_id"],
                "updated_widgets": msg["updated_widgets"],
            }

        elif msg["type"] == "output_value" and "cell_id" in msg:
            await gql_query(
                auth=auth,
                query="""
                    mutation UpdateCellOutputValue($id: BigInt!, $data: String!) {
                        updatePlotTransformInfo(
                            input: { id: $id, patch: { cellViewerData: $data } }
                        ) {
                            clientMutationId
                        }
                    }
                """,
                variables={"id": msg["cell_id"], "data": orjson.dumps(msg).decode()},
            )

            msg = {"type": msg["type"], "cell_id": msg["cell_id"]}

        elif msg["type"] == "output_value" and "viewer_id" in msg:
            await gql_query(
                auth=auth,
                query="""
                    mutation UpdatePlotCellValueViewer($id: BigInt!, $data: String!) {
                        updatePlotCellValueViewer(
                            input: { id: $id, patch: { cellViewerData: $data } }
                        ) {
                            clientMutationId
                        }
                    }
                """,
                variables={"id": msg["viewer_id"], "data": orjson.dumps(msg).decode()},
            )

            msg = {"type": msg["type"], "viewer_id": msg["viewer_id"]}

        elif msg["type"] == "plot_data":
            await gql_query(
                auth=auth,
                query="""
                    mutation UpdatePlotInfo($id: BigInt!, $data: String!) {
                        updatePlotInfo(
                            input: { id: $id, patch: { transformData: $data } }
                        ) {
                            clientMutationId
                        }
                    }
                """,
                variables={"id": msg["plot_id"], "data": orjson.dumps(msg).decode()},
            )

            msg = {"type": msg["type"], "plot_id": msg["plot_id"]}

        await broadcast_message(orjson.dumps(msg).decode())


async def handle_kernel_io(stream: asyncio.StreamReader, *, name: str) -> None:
    await ready_ev.wait()
    print(f"Strating kernel {name} listener")
    while True:
        data = await stream.read(4096)
        print(f"[kernel] ({stream}): {data}")
        if len(data) == 0:
            break

        # todo(rteqs): message are getting sent out of order. kernel stdio is sent then start_cell
        # which clears the initial output / logs. We need to fully rework logging, as we have another
        # hack in the kernel related to logs.
        await asyncio.sleep(0.1)
        await broadcast_message(
            orjson.dumps(
                {
                    "type": "kernel_stdio",
                    "active_cell": active_cell,
                    "stream": name,
                    "data": data.decode(errors="replace"),
                }
            ).decode()
        )


async def start_kernel_proc() -> None:
    await add_pod_event(auth=auth_token_sdk, event_type="runtime_starting")
    conn_k = k_proc.conn_k = await SocketIo.from_socket(sock)
    async_tasks.append(
        asyncio.create_task(handle_kernel_messages(k_proc.conn_k, auth_token_sdk))
    )

    print("Starting kernel subprocess")
    proc = k_proc.proc = await asyncio.create_subprocess_exec(
        sys.executable,
        (dir_p / "kernel.py"),
        str(sock_k_fd),
        pass_fds=[sock_k_fd],
        stdin=asyncio.subprocess.DEVNULL,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        preexec_fn=lambda: os.nice(1),
    )
    assert proc.stdout is not None
    assert proc.stderr is not None

    async_tasks.extend(
        [
            asyncio.create_task(handle_kernel_io(proc.stdout, name="stdout")),
            asyncio.create_task(handle_kernel_io(proc.stderr, name="stderr")),
        ]
    )

    resp = await gql_query(
        auth=auth_token_sdk,
        query="""
            query plotsNotebookKernelState($pod_id: BigInt!) {
                plotsNotebookKernelState(argPodId: $pod_id)
            }
        """,
        variables={"pod_id": pod_id},
    )

    data = validate(resp, PlotsNotebookKernelStateResp)
    k_state = data.data.plotsNotebookKernelState

    await add_pod_event(auth=auth_token_sdk, event_type="runtime_ready")
    await ready_ev.wait()
    # todo(rteqs): need to init plot configs
    await conn_k.send(
        {
            "type": "init",
            "widget_states": k_state.widget_states,
            "cell_output_selections": k_state.cell_output_selections,
            "plot_data_selections": k_state.plot_data_selections,
            "viewer_cell_data": k_state.viewer_cell_data,
            "plot_configs": k_state.plot_configs,
        }
    )


async def stop_kernel_proc() -> None:
    ready_ev.clear()

    proc = k_proc.proc
    if proc is not None:
        try:
            proc.terminate()
            await asyncio.wait_for(proc.wait(), timeout=10)
        except TimeoutError:
            print("Error terminating kernel process")
            proc.kill()
            await proc.wait()

    for task in async_tasks:
        task.cancel()


async def shutdown() -> None:
    await stop_kernel_proc()
    with contextlib.suppress(Exception):
        sock_k.close()

    with contextlib.suppress(Exception):
        sock.close()

    with contextlib.suppress(Exception):
        sess = get_global_http_sess()
        await sess.close()
