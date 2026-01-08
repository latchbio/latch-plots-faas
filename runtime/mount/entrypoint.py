import asyncio
import contextlib
import os
import socket
import sys
import traceback
from asyncio.subprocess import Process
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import IO, TypedDict, TypeVar

import orjson
from latch_asgi.context.websocket import Context
from latch_data_validation.data_validation import validate

from runtime.mount.plots_context_manager import PlotsContextManager

from .headless_browser import HeadlessBrowser
from .socketio import SocketIo
from .utils import (
    KernelSnapshotStatus,
    PlotConfig,
    get_global_http_sess,
    gql_query,
    orjson_encoder,
)

dir_p = Path(__file__).parent


T = TypeVar("T")


sock, sock_k = socket.socketpair(family=socket.AF_UNIX)
sock.setblocking(False)
sock_k_fd = sock_k.detach()

sock_a, sock_agent = socket.socketpair(family=socket.AF_UNIX)
sock_a.setblocking(False)
sock_agent_fd = sock_agent.detach()

ready_ev = asyncio.Event()

active_cell: str | None = None


class CellOutputs(TypedDict):
    outputs: list[str]
    dataframe_outputs: list[str]
    figure_outputs: list[str]
    static_figure_outputs: list[str]


cell_status: dict[str, str] = {}
cell_sequencers: dict[str, int] = {}
cell_last_run_outputs: dict[str, CellOutputs] = {}


@dataclass
class KernelSnapshotState:
    status: KernelSnapshotStatus = "done"


kernel_snapshot_state = KernelSnapshotState()

async_tasks: list[asyncio.Task] = []


latch_p = Path("/root/.latch")
sdk_token = (latch_p / "token").read_text()
auth_token_sdk = f"Latch-SDK-Token {sdk_token}"

pod_id = int((latch_p / "id").read_text())
pod_session_id = (latch_p / "session-id").read_text()

plots_ctx_manager = PlotsContextManager()
current_agent_ctx: Context | None = None


@dataclass
class KernelProc:
    conn_k: SocketIo | None = None
    proc: Process | None = None


k_proc = KernelProc()


@dataclass
class AgentProc:
    conn_a: SocketIo | None = None
    proc: Process | None = None
    log_file: IO[str] | None = None


a_proc = AgentProc()

headless_browser: HeadlessBrowser | None = None
latest_local_storage: dict[str, str] | None = None


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


@dataclass(frozen=True)
class PlotUpsertValueViewerCreateResp:
    bigInt: str | None


@dataclass(frozen=True)
class PlotUpsertValueViewerDataResp:
    upsertPlotCellValueViewer: PlotUpsertValueViewerCreateResp


@dataclass(frozen=True)
class PlotUpsertValueViewerGQLResp:
    data: PlotUpsertValueViewerDataResp


@dataclass(frozen=True)
class TmpPlotsNotebookKernelSnapshotMode:
    tmpPlotsNotebookKernelSnapshotMode: bool


@dataclass(frozen=True)
class TmpPlotsNotebookKernelSnapshotModeResp:
    data: TmpPlotsNotebookKernelSnapshotMode


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
    except Exception as e:
        if (
            'duplicate key value violates unique constraint "pod_session_events_unique_sess_id_event_type"'
            in str(e)
        ):
            # todo(maximsmol): fix this properly
            pass

        traceback.print_exc()


async def handle_kernel_messages(conn_k: SocketIo, auth: str) -> None:
    global active_cell

    print("Starting kernel message listener")
    while True:
        msg = await conn_k.recv()

        try:
            # note(aidan): h5 messages can be mbs and stuffing into journal fills the disc
            if msg["type"] != "h5":
                print(">", msg)

            if msg["type"] == "ready":
                ready_ev.set()
                await add_pod_event(auth=auth, event_type="kernel_ready")
                continue

            if msg["type"] == "debug_state":
                msg["sess_hash"] = pod_session_id

            elif msg["type"] == "start_cell":
                cell_id = msg["cell_id"]

                if cell_id is not None:
                    # note: active_cell is set to None when on_tick_finished is called to tell console to mark all cells as ran.
                    # active_cell is not overwritten here becuase kernel_stdio messages require an active_cell field,
                    # if active_cell is None in kernel_stdio, console will display an error
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
                outputs_data: CellOutputs = {
                    "outputs": msg["outputs"],
                    "dataframe_outputs": msg["dataframe_outputs"],
                    "figure_outputs": msg["figure_outputs"],
                    "static_figure_outputs": msg["static_figure_outputs"],
                }
                cell_last_run_outputs[cell_id] = outputs_data

                exc = msg.get("exception")
                if exc is not None:
                    exc = orjson.dumps({"string": exc[-9000:]}).decode()

                await gql_query(
                    auth=auth,
                    query="""
                        mutation UpdateCellResult($id: BigInt!, $exception: String,  $data: String) {
                            updatePlotTransformInfo(
                                input: { id: $id, patch: { exception: $exception, outputsData: $data } }
                            ) {
                                clientMutationId
                            }
                        }
                    """,
                    variables={
                        "id": msg["cell_id"],
                        "exception": exc,
                        "data": orjson.dumps(outputs_data).decode(),
                    },
                )

                msg = {
                    "type": msg["type"],
                    "cell_id": msg["cell_id"],
                    "has_exception": exc is not None,
                    "exception": exc,
                    **outputs_data,
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
                    variables={
                        "id": msg["cell_id"],
                        "data": orjson.dumps(
                            msg,
                            option=orjson.OPT_SERIALIZE_NUMPY,
                            default=orjson_encoder,
                        ).decode(),
                    },
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
                    variables={
                        "id": msg["viewer_id"],
                        "data": orjson.dumps(
                            msg,
                            option=orjson.OPT_SERIALIZE_NUMPY,
                            default=orjson_encoder,
                        ).decode(),
                    },
                )

                msg = {"type": msg["type"], "viewer_id": msg["viewer_id"]}

            elif msg["type"] == "cell_value_viewer_init":
                assert plots_ctx_manager.notebook_id is not None

                resp = await gql_query(
                    auth=auth,
                    query="""
                        mutation UpsertPlotCellValueViewer(
                            $notebookId: BigInt!,
                            $widgetConnectionKey: String!
                        ) {
                            upsertPlotCellValueViewer(
                                input: {
                                    argNotebookId: $notebookId,
                                    argWidgetConnectionKey: $widgetConnectionKey
                                }
                            ) {
                                bigInt
                            }
                        }
                    """,
                    variables={
                        "notebookId": plots_ctx_manager.notebook_id,
                        "widgetConnectionKey": msg["value_viewer_key"],
                    },
                )

                data = validate(resp, PlotUpsertValueViewerGQLResp)
                viewer_id = data.data.upsertPlotCellValueViewer.bigInt
                if viewer_id is not None:
                    await conn_k.send(
                        {
                            "type": "get_global",
                            "viewer_id": viewer_id,
                            "key": msg["global_key"],
                        }
                    )

                continue

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
                    variables={
                        "id": msg["plot_id"],
                        "data": orjson.dumps(
                            msg,
                            option=orjson.OPT_SERIALIZE_NUMPY,
                            default=orjson_encoder,
                        ).decode(),
                    },
                )

                msg = {"type": msg["type"], "plot_id": msg["plot_id"]}

            elif msg["type"] in {"load_kernel_snapshot", "save_kernel_snapshot"}:
                kernel_snapshot_state.status = msg["status"]

            elif msg["type"] == "globals_summary" and "agent_tx_id" in msg:
                tx_id = msg.get("agent_tx_id")

                if a_proc.conn_a is not None:
                    print(f"[entrypoint] Routing globals response to agent (tx_id={tx_id})")
                    await a_proc.conn_a.send({
                        "type": "agent_action_response",
                        "tx_id": tx_id,
                        "status": "success",
                        "summary": msg.get("summary", {})
                    })
                else:
                    print("[entrypoint] Could not route globals response: agent not connected")

                continue

            elif msg["type"] == "reactivity_summary" and "agent_tx_id" in msg:
                tx_id = msg.get("agent_tx_id")

                if a_proc.conn_a is not None:
                    print(f"[entrypoint] Routing reactivity response to agent (tx_id={tx_id})")
                    await a_proc.conn_a.send({
                        "type": "agent_action_response",
                        "tx_id": tx_id,
                        "status": "success",
                        "cell_reactivity": msg.get("cell_reactivity", {}),
                    })
                    continue

            elif msg["type"] == "execute_code_response" and "agent_tx_id" in msg:
                tx_id = msg.get("agent_tx_id")

                if a_proc.conn_a is not None:
                    print(f"[entrypoint] Routing execute_code response to agent (tx_id={tx_id})")
                    await a_proc.conn_a.send({
                        "type": "agent_action_response",
                        "tx_id": tx_id,
                        "status": msg.get("status", "error"),
                        "result": msg.get("result"),
                        "stdout": msg.get("stdout"),
                        "stderr": msg.get("stderr"),
                        "exception": msg.get("exception"),
                        "error": msg.get("error")
                    })
                    continue

            elif msg["type"] == "get_global_info_response" and "agent_tx_id" in msg:
                tx_id = msg.get("agent_tx_id")

                if a_proc.conn_a is not None:
                    print(f"[entrypoint] Routing get_global_info response to agent (tx_id={tx_id})")
                    await a_proc.conn_a.send({
                        "type": "agent_action_response",
                        "tx_id": tx_id,
                        "status": msg.get("status", "error"),
                        "info": msg.get("info"),
                        "error": msg.get("error")
                    })
                    continue

            await plots_ctx_manager.broadcast_message(orjson.dumps(msg).decode())

        except Exception:
            err_msg = {"type": "error", "data": traceback.format_exc()}
            await plots_ctx_manager.broadcast_message(orjson.dumps(err_msg).decode())


async def handle_agent_messages(conn_a: SocketIo) -> None:
    ts = lambda: datetime.utcnow().isoformat() + "Z"
    print(f"{ts()} [entrypoint] Starting agent message listener")
    while True:
        msg = await conn_a.recv()
        msg_type = msg.get("type", "unknown")
        tx_id = msg.get("tx_id")
        action = msg.get("action")

        if msg_type != "agent_stream_delta":
            print(f"{ts()} [entrypoint] Agent > {msg_type}")

        if msg_type == "agent_action" and msg.get("action") == "request_reactivity_summary":
            if k_proc.conn_k is not None:
                print(f"{ts()} [entrypoint] Routing reactivity request to kernel (tx_id={tx_id})")

                await k_proc.conn_k.send({
                    "type": "reactivity_summary",
                    "agent_tx_id": tx_id
                })
            else:
                print("[entrypoint] Kernel not connected, cannot route reactivity request")
                await conn_a.send({
                    "type": "agent_action_response",
                    "tx_id": msg.get("tx_id"),
                    "status": "error",
                    "error": "Kernel not connected"
                })
            continue

        if msg_type == "agent_action" and msg.get("action") == "execute_code":
            if k_proc.conn_k is not None:
                code = msg.get("params", {}).get("code", "")
                print(f"{ts()} [entrypoint] Routing execute_code to kernel (tx_id={tx_id})")

                await k_proc.conn_k.send({
                    "type": "execute_code",
                    "code": code,
                    "agent_tx_id": tx_id
                })
            else:
                await conn_a.send({
                    "type": "agent_action_response",
                    "tx_id": msg.get("tx_id"),
                    "status": "error",
                    "error": "Kernel not connected"
                })
            continue

        if msg_type == "agent_action" and msg.get("action") == "get_global_info":
            if k_proc.conn_k is not None:
                key = msg.get("params", {}).get("key", "")
                print(f"{ts()} [entrypoint] Routing get_global_info to kernel (tx_id={tx_id}, key={key})")

                await k_proc.conn_k.send({
                    "type": "get_global_info",
                    "key": key,
                    "agent_tx_id": tx_id
                })
            else:
                await conn_a.send({
                    "type": "agent_action_response",
                    "tx_id": msg.get("tx_id"),
                    "status": "error",
                    "error": "Kernel not connected"
                })
            continue

        if current_agent_ctx is None:
            print(f"{ts()} [entrypoint] No websocket client connected, skipping message: {msg_type} action={action} tx_id={tx_id}")
            continue

        try:
            print(f"{ts()} [entrypoint] Forwarding agent message to browser ctx action={action} tx_id={tx_id}")
            await current_agent_ctx.send_message(orjson.dumps(msg).decode())
        except Exception as e:
            print(f"{ts()} [entrypoint] Error forwarding message to browser ctx action={action} tx_id={tx_id}: {e}")


async def start_kernel_proc() -> None:
    await add_pod_event(auth=auth_token_sdk, event_type="runtime_starting")
    conn_k = k_proc.conn_k = await SocketIo.from_socket(sock)
    async_tasks.append(
        asyncio.create_task(handle_kernel_messages(k_proc.conn_k, auth_token_sdk))
    )

    print("Starting kernel subprocess")
    k_proc.proc = await asyncio.create_subprocess_exec(
        sys.executable,
        (dir_p / "kernel.py"),
        str(sock_k_fd),
        pass_fds=[sock_k_fd],
        stdin=asyncio.subprocess.DEVNULL,
        preexec_fn=lambda: os.nice(1),
    )

    k_state: KernelState | None = None
    try:
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
    except Exception:
        err_msg = {"type": "error", "data": traceback.format_exc()}
        await plots_ctx_manager.broadcast_message(orjson.dumps(err_msg).decode())
        traceback.print_exc()

    if k_state is None:
        return

    # todo(kenny): separate query for backwards compatability. Pull into main
    # fn when "snapshot mode" merged and works well
    try:
        resp = await gql_query(
            auth=auth_token_sdk,
            query="""
                query tmpPlotsNotebookKernelSnapshotMode($pod_id: BigInt!) {
                    tmpPlotsNotebookKernelSnapshotMode(argPodId: $pod_id)
                }
            """,
            variables={"pod_id": pod_id},
        )

        data = validate(resp, TmpPlotsNotebookKernelSnapshotModeResp)
        session_snapshot_mode = data.data.tmpPlotsNotebookKernelSnapshotMode
        if session_snapshot_mode:
            kernel_snapshot_state.status = "loading"
    except Exception:
        err_msg = {"type": "error", "data": traceback.format_exc()}
        await plots_ctx_manager.broadcast_message(orjson.dumps(err_msg).decode())
        session_snapshot_mode = False

    await add_pod_event(auth=auth_token_sdk, event_type="runtime_ready")
    await ready_ev.wait()
    await conn_k.send(
        {
            "type": "init",
            "widget_states": k_state.widget_states,
            "cell_output_selections": k_state.cell_output_selections,
            "plot_data_selections": k_state.plot_data_selections,
            "viewer_cell_data": k_state.viewer_cell_data,
            "plot_configs": k_state.plot_configs,
            "session_snapshot_mode": session_snapshot_mode,
        }
    )


async def start_agent_proc() -> None:
    conn_a = a_proc.conn_a = await SocketIo.from_socket(sock_a)
    async_tasks.append(
        asyncio.create_task(handle_agent_messages(a_proc.conn_a))
    )

    log_path = Path("/var/log/agent.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    a_proc.log_file = open(log_path, "a")

    print("Starting agent subprocess")
    a_proc.proc = await asyncio.create_subprocess_exec(
        sys.executable,
        (dir_p / "agent.py"),
        str(sock_agent_fd),
        pass_fds=[sock_agent_fd],
        stdin=asyncio.subprocess.DEVNULL,
        stdout=a_proc.log_file,
        stderr=a_proc.log_file,
        preexec_fn=lambda: os.nice(5),
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


async def stop_agent_proc() -> None:
    proc = a_proc.proc
    if proc is not None:
        if proc.returncode is None:
            try:
                proc.terminate()
                await asyncio.wait_for(proc.wait(), timeout=2)
            except TimeoutError:
                print("Error terminating agent process")
                proc.kill()
                await proc.wait()
            except ProcessLookupError:
                pass
        else:
            print(f"[entrypoint] Agent process already exited with code {proc.returncode}")

    if a_proc.log_file is not None:
        a_proc.log_file.close()


async def start_headless_browser(notebook_id: str, local_storage: dict[str, str]) -> None:
    global headless_browser

    if headless_browser is not None:
        print("[entrypoint] Headless browser already running")
        return

    try:
        notebook_url = f"https://console.latch.bio/plots/{notebook_id}"

        print(f"[entrypoint] Starting headless browser for notebook {notebook_id}")
        headless_browser = HeadlessBrowser()
        await headless_browser.start(
            notebook_url,
            local_storage=local_storage,
        )
        with contextlib.suppress(Exception):
            await headless_browser.screenshot("/var/log/plots-headless.png")
            print("[entrypoint] Saved headless browser screenshot to /var/log/plots-headless.png")
        print(f"[entrypoint] Headless browser ready for notebook {notebook_id}")

    except Exception as e:
        print(f"[entrypoint] Error starting headless browser: {e}")
        traceback.print_exc()


async def shutdown() -> None:
    await stop_kernel_proc()
    await stop_agent_proc()

    global headless_browser
    if headless_browser is not None:
        await headless_browser.stop()
        headless_browser = None

    with contextlib.suppress(Exception):
        sock_k.close()
    with contextlib.suppress(Exception):
        sock_agent.close()

    with contextlib.suppress(Exception):
        sock.close()
    with contextlib.suppress(Exception):
        sock_a.close()

    with contextlib.suppress(Exception):
        sess = get_global_http_sess()
        await sess.close()
