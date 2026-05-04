import asyncio
import contextlib
import json
import os
import socket
import sys
import traceback
from asyncio.subprocess import Process
from collections.abc import Awaitable, Callable
from dataclasses import dataclass, field
from io import TextIOWrapper
from pathlib import Path
from typing import TypedDict, TypeVar

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

ready_ev = asyncio.Event()


class CellOutputs(TypedDict):
    outputs: list[str]
    dataframe_outputs: list[str]
    figure_outputs: list[str]
    static_figure_outputs: list[str]


cell_status: dict[str, str] = {}
cell_sequencers: dict[str, int] = {}
cell_last_run_outputs: dict[str, CellOutputs] = {}
latest_reactivity_summary: dict[str, dict[str, list[str]]] | None = None


@dataclass
class KernelSnapshotState:
    status: KernelSnapshotStatus = "done"


kernel_snapshot_state = KernelSnapshotState()

async_tasks: list[asyncio.Task[object]] = []


latch_p = Path("/root/.latch")
sdk_token = (latch_p / "token").read_text()
auth_token_sdk = f"Latch-SDK-Token {sdk_token}"

pod_id = int((latch_p / "id").read_text())
pod_session_id = (latch_p / "session-id").read_text()

# notebook_id is directly provided in the case of the eval harness (no corresponding pod)
notebook_id_file = latch_p / "notebook-id"
notebook_id = None
with contextlib.suppress(Exception):
    notebook_id = (
        notebook_id_file.read_text().strip() if notebook_id_file.exists() else None
    )

plots_ctx_manager = PlotsContextManager()

user_agent_ctx: Context | None = None

# backend browser
action_handler_ctx: Context | None = None
action_handler_ready_ev = asyncio.Event()

pending_user_browser_actions: dict[str, dict] = {}  # tx_id -> original message


def mark_action_handled(tx_id: str) -> None:
    _ = pending_user_browser_actions.pop(tx_id, None)


async def handle_user_disconnect_fallback() -> None:
    actions_to_fallback = list(pending_user_browser_actions.items())
    pending_user_browser_actions.clear()

    for tx_id, msg in actions_to_fallback:
        if action_handler_ctx is not None:
            await action_handler_ctx.send_message(orjson.dumps(msg).decode())
        elif a_proc.msg_io is not None:
            await a_proc.msg_io.send({
                "type": "agent_action_response",
                "tx_id": tx_id,
                "status": "error",
                "error": "User disconnected and backend browser not available",
            })


shutting_down = False


@dataclass(kw_only=True)
class ProcessManager:
    name: str
    proc: Process | None = None
    msg_io: SocketIo | None = None
    log_io: TextIOWrapper | None = None

    started: asyncio.Condition = field(default_factory=asyncio.Condition)

    async def create_msg_io(self) -> int:
        sock, sock_k = socket.socketpair(family=socket.AF_UNIX)
        sock.setblocking(False)  # noqa: FBT003
        sock_k_fd = sock_k.detach()

        self.msg_io = await SocketIo.from_socket(sock)
        return sock_k_fd

    async def wait_for_start(self) -> Process:
        async with self.started:
            while True:
                if self.proc is not None and self.proc.returncode is None:
                    return self.proc

                _ = await self.started.wait()

    async def stop(self) -> None:
        if self.proc is not None and self.proc.returncode is None:
            try:
                try:
                    self.proc.terminate()
                    _ = await asyncio.wait_for(self.proc.wait(), timeout=2)
                except TimeoutError:
                    print(f"Process did not terminate in 2 seconds: {self.name}")
                    self.proc.kill()
                    _ = await self.proc.wait()
                except ProcessLookupError:
                    pass
            except Exception:
                traceback.print_exc()

        if self.log_io is not None:
            try:
                self.log_io.close()
            except Exception:
                traceback.print_exc()
        self.log_io = None

        if self.msg_io is not None:
            try:
                self.msg_io.sock.shutdown(socket.SHUT_RDWR)
                self.msg_io.sock.close()
            except Exception:
                traceback.print_exc()
        self.msg_io = None

    async def watchdog(
        self, *, cleanup: Callable[..., Awaitable[object]] | None = None
    ) -> None:
        while True:
            print(f"[watchdog] <{self.name}> Starting")
            proc = await self.wait_for_start()

            code = await proc.wait()
            if shutting_down:
                return

            print(f"[watchdog] <{self.name}> Exited with code {code}")
            await plots_ctx_manager.broadcast_message(
                orjson.dumps({
                    "type": "managed_process_exit",
                    "name": self.name,
                    "exit_code": code,
                }).decode()
            )

            if cleanup is not None:
                _ = await cleanup()


k_proc = ProcessManager(name="kernel")
a_proc = ProcessManager(name="agent")

headless_browser: HeadlessBrowser | None = None
latest_local_storage: dict[str, str] | None = None
headless_browser_notebook_id: str | None = None


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
class PlotsNotebookKernelStateByNotebook:
    plotsNotebookKernelStateByNotebook: KernelState


@dataclass(frozen=True)
class PlotsNotebookKernelStateByNotebookResp:
    data: PlotsNotebookKernelStateByNotebook


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


@dataclass(frozen=True)
class AgentBootstrapPodInfoNode:
    plotNotebookId: str | None
    accountId: str | None


@dataclass(frozen=True)
class AgentBootstrapPodInfos:
    nodes: list[AgentBootstrapPodInfoNode]


@dataclass(frozen=True)
class AgentBootstrapPodInfoData:
    podInfos: AgentBootstrapPodInfos


@dataclass(frozen=True)
class AgentBootstrapPodInfoResp:
    data: AgentBootstrapPodInfoData


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
    global latest_reactivity_summary

    print("Starting kernel message listener")
    while True:
        msg_raw = await conn_k.recv_raw()
        msg = orjson.loads(msg_raw)

        try:
            print(f"> {msg.get('type', 'unknown')}: {len(msg_raw)}")

            if msg["type"] == "ready":
                ready_ev.set()
                await add_pod_event(auth=auth, event_type="kernel_ready")
                continue

            if msg["type"] == "debug_state":
                msg["sess_hash"] = pod_session_id

            elif msg["type"] == "start_cell":
                cell_id = msg["cell_id"]

                if cell_id is not None:
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

                data = msg.get("source")
                if data is not None:
                    data = orjson.dumps(
                        data, option=orjson.OPT_SERIALIZE_NUMPY, default=orjson_encoder
                    ).decode()

                resp = await gql_query(
                    auth=auth,
                    query="""
                        mutation UpsertPlotCellValueViewer(
                            $notebookId: BigInt!,
                            $widgetConnectionKey: String!,
                            $data: String
                        ) {
                            upsertPlotCellValueViewer(
                                input: {
                                    argNotebookId: $notebookId,
                                    argWidgetConnectionKey: $widgetConnectionKey,
                                    argViewerData: $data
                                }
                            ) {
                                bigInt
                            }
                        }
                    """,
                    variables={
                        "notebookId": plots_ctx_manager.notebook_id,
                        "widgetConnectionKey": msg["value_viewer_key"],
                        "data": data,
                    },
                )

                data = validate(resp, PlotUpsertValueViewerGQLResp)
                viewer_id = data.data.upsertPlotCellValueViewer.bigInt

                if viewer_id is None:
                    continue

                if msg["global_key"] is not None:
                    await conn_k.send({
                        "type": "get_global",
                        "viewer_id": viewer_id,
                        "key": msg["global_key"],
                    })

                    continue

                msg = {"type": "output_value", "viewer_id": viewer_id}

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

            elif msg["type"] == "reactivity_summary_update":
                cell_reactivity = msg.get("cell_reactivity")
                if isinstance(cell_reactivity, dict):
                    latest_reactivity_summary = cell_reactivity
                continue

            elif msg["type"] == "globals_summary" and "agent_tx_id" in msg:
                tx_id = msg.get("agent_tx_id")

                if a_proc.msg_io is not None:
                    print(
                        f"[entrypoint] Routing globals response to agent (tx_id={tx_id})"
                    )
                    await a_proc.msg_io.send({
                        "type": "agent_action_response",
                        "tx_id": tx_id,
                        "status": "success",
                        "summary": msg.get("summary", {}),
                    })
                else:
                    print(
                        "[entrypoint] Could not route globals response: agent not connected"
                    )

                continue

            elif msg["type"] == "execute_code_response" and "agent_tx_id" in msg:
                tx_id = msg.get("agent_tx_id")

                if a_proc.msg_io is not None:
                    print(
                        f"[entrypoint] Routing execute_code response to agent (tx_id={tx_id})"
                    )
                    await a_proc.msg_io.send({
                        "type": "agent_action_response",
                        "tx_id": tx_id,
                        "status": msg.get("status", "error"),
                        "result": msg.get("result"),
                        "stdout": msg.get("stdout"),
                        "stderr": msg.get("stderr"),
                        "exception": msg.get("exception"),
                        "error": msg.get("error"),
                    })
                    continue

            elif msg["type"] == "get_global_info_response" and "agent_tx_id" in msg:
                tx_id = msg.get("agent_tx_id")

                if a_proc.msg_io is not None:
                    print(
                        f"[entrypoint] Routing get_global_info response to agent (tx_id={tx_id})"
                    )
                    await a_proc.msg_io.send({
                        "type": "agent_action_response",
                        "tx_id": tx_id,
                        "status": msg.get("status", "error"),
                        "info": msg.get("info"),
                        "error": msg.get("error"),
                    })
                    continue

            await plots_ctx_manager.broadcast_message(orjson.dumps(msg).decode())

        except Exception:
            err_msg = {"type": "error", "data": traceback.format_exc()}
            await plots_ctx_manager.broadcast_message(orjson.dumps(err_msg).decode())


async def handle_agent_messages(conn_a: SocketIo) -> None:
    print("[entrypoint] Starting agent message listener")
    while True:
        msg = await conn_a.recv()
        msg_type = msg.get("type", "unknown")
        tx_id = msg.get("tx_id")
        action = msg.get("action")

        if msg_type != "agent_stream_delta" and "agent_stream_block" not in msg_type:
            print(f"[entrypoint] Agent > {msg_type}")

        if msg_type == "agent_action":
            print(f"[entrypoint]       > (action={action})")

        if (
            msg_type == "agent_action"
            and msg.get("action") == "request_reactivity_summary"
        ):
            if latest_reactivity_summary is None:
                print("[entrypoint] Reactivity summary not ready yet")
                await conn_a.send({
                    "type": "agent_action_response",
                    "tx_id": msg.get("tx_id"),
                    "status": "error",
                    "error": "Reactivity summary not ready yet",
                })
            else:
                await conn_a.send({
                    "type": "agent_action_response",
                    "tx_id": msg.get("tx_id"),
                    "status": "success",
                    "cell_reactivity": latest_reactivity_summary,
                })
            continue

        if msg_type == "agent_action" and msg.get("action") == "execute_code":
            if k_proc.msg_io is not None:
                code = msg.get("params", {}).get("code", "")

                await k_proc.msg_io.send({
                    "type": "execute_code",
                    "code": code,
                    "agent_tx_id": tx_id,
                })
            else:
                await conn_a.send({
                    "type": "agent_action_response",
                    "tx_id": msg.get("tx_id"),
                    "status": "error",
                    "error": "Kernel not connected",
                })
            continue

        if msg_type == "agent_action" and msg.get("action") == "get_global_info":
            if k_proc.msg_io is not None:
                key = msg.get("params", {}).get("key", "")

                await k_proc.msg_io.send({
                    "type": "get_global_info",
                    "key": key,
                    "agent_tx_id": tx_id,
                })
            else:
                await conn_a.send({
                    "type": "agent_action_response",
                    "tx_id": msg.get("tx_id"),
                    "status": "error",
                    "error": "Kernel not connected",
                })
            continue

        if msg_type == "agent_action":
            force_backend_browser_retry = msg.get("params", {}).get(
                "force_backend_browser_retry", False
            )

            if action in {
                "smart_ui_spotlight",
                "h5_filter_by",
                "h5_color_by",
                "h5_set_selected_obsm_key",
                "h5_set_background_image",
                "h5_open_image_aligner",
                "h5_autoscale",
                "h5_zoom",
                "h5_set_background_image_visibility",
                "h5_add_selected_cells_to_categorical_obs",
                "h5_set_marker_opacity",
                "h5_manage_obs",
            }:
                if user_agent_ctx is not None:
                    target_ctx = user_agent_ctx
                    track_for_disconnect = False
                else:
                    print(
                        f"[entrypoint] User browser unavailable for user-specific action: {action}"
                    )
                    await conn_a.send({
                        "type": "agent_action_response",
                        "tx_id": tx_id,
                        "status": "error",
                        "error": (
                            """
                            The agent cannot interact with the user as they are not actively viewing
                            the page but this action requires user interaction. Please proceed
                            without relying on user input, or inform the user they need to be on the
                            page to proceed.
                            """
                        ),
                    })
                    continue
            elif force_backend_browser_retry:
                if action_handler_ctx is not None:
                    target_ctx = action_handler_ctx
                    track_for_disconnect = False
                else:
                    print(
                        f"[entrypoint] Backend browser not connected for forced action: {action}"
                    )
                    await conn_a.send({
                        "type": "agent_action_response",
                        "tx_id": tx_id,
                        "status": "error",
                        "error": "Backend browser not connected. Cannot execute browser actions.",
                    })
                    continue
            elif user_agent_ctx is not None:
                target_ctx = user_agent_ctx
                track_for_disconnect = True
            elif action_handler_ctx is not None:
                target_ctx = action_handler_ctx
                track_for_disconnect = False
            else:
                print(f"[entrypoint] No browser connected for agent_action: {action}")
                await conn_a.send({
                    "type": "agent_action_response",
                    "tx_id": tx_id,
                    "status": "error",
                    "error": "No browser connected. Cannot execute browser actions.",
                })
                continue
        else:
            target_ctx = user_agent_ctx
            track_for_disconnect = False
            if target_ctx is None:
                target_ctx = action_handler_ctx
            if target_ctx is None:
                continue

        try:
            if track_for_disconnect and tx_id is not None:
                pending_user_browser_actions[tx_id] = msg

            await target_ctx.send_message(orjson.dumps(msg).decode())
        except Exception as e:
            print(f"[entrypoint] Error forwarding {msg_type} message: {e}")
            if track_for_disconnect and tx_id is not None:
                mark_action_handled(tx_id)
                if action_handler_ctx is not None:
                    await action_handler_ctx.send_message(orjson.dumps(msg).decode())
                elif a_proc.msg_io is not None:
                    await a_proc.msg_io.send({
                        "type": "agent_action_response",
                        "tx_id": tx_id,
                        "status": "error",
                        "error": f"Failed to forward action: {e}",
                    })


async def start_kernel_proc() -> None:
    ready_ev.clear()

    cell_status.clear()
    cell_sequencers.clear()
    cell_last_run_outputs.clear()

    global latest_reactivity_summary

    latest_reactivity_summary = None

    kernel_snapshot_state.status = "done"

    await add_pod_event(auth=auth_token_sdk, event_type="runtime_starting")

    sock_k_fd = await k_proc.create_msg_io()
    assert k_proc.msg_io is not None

    async_tasks.append(
        asyncio.create_task(handle_kernel_messages(k_proc.msg_io, auth_token_sdk))
    )

    # todo(maximsmol): log rotation? syslog?
    log_path = Path("/var/log/kernel.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    a_proc.log_io = Path(log_path).open("a", encoding="utf-8")

    print("Starting kernel subprocess")
    k_proc.proc = await asyncio.create_subprocess_exec(
        sys.executable,
        (dir_p / "kernel.py"),
        str(sock_k_fd),
        pass_fds=[sock_k_fd],
        stdin=asyncio.subprocess.DEVNULL,
        stdout=a_proc.log_io,
        stderr=a_proc.log_io,
        preexec_fn=lambda: os.nice(1),
    )
    async with k_proc.started:
        k_proc.started.notify_all()

    _ = Path(f"/proc/{k_proc.proc.pid}/oom_score_adj").write_text(
        "100\n", encoding="utf-8"
    )

    k_state: KernelState | None = None
    try:
        # eval harness case
        if notebook_id is not None:
            resp = await gql_query(
                auth=auth_token_sdk,
                query="""
                    query plotsNotebookKernelStateByNotebook($notebook_id: BigInt!) {
                        plotsNotebookKernelStateByNotebook(argNotebookId: $notebook_id)
                    }
                """,
                variables={"notebook_id": notebook_id},
            )
            data = validate(resp, PlotsNotebookKernelStateByNotebookResp)
            k_state = data.data.plotsNotebookKernelStateByNotebook
        else:
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
            kernel_snapshot_state.status = "start"
    except Exception:
        err_msg = {"type": "error", "data": traceback.format_exc()}
        await plots_ctx_manager.broadcast_message(orjson.dumps(err_msg).decode())
        session_snapshot_mode = False

    await add_pod_event(auth=auth_token_sdk, event_type="runtime_ready")
    _ = await ready_ev.wait()
    await k_proc.msg_io.send({
        "type": "init",
        "widget_states": k_state.widget_states,
        "cell_output_selections": k_state.cell_output_selections,
        "plot_data_selections": k_state.plot_data_selections,
        "viewer_cell_data": k_state.viewer_cell_data,
        "plot_configs": k_state.plot_configs,
        "session_snapshot_mode": session_snapshot_mode,
        "kernel_pid": k_proc.proc.pid,
    })


async def start_agent_proc() -> None:
    # todo(maximsmol): unify with start_kernel_proc

    sock_agent_fd = await a_proc.create_msg_io()
    assert a_proc.msg_io is not None

    async_tasks.append(asyncio.create_task(handle_agent_messages(a_proc.msg_io)))

    log_path = Path("/var/log/agent.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    a_proc.log_io = Path(log_path).open("a", encoding="utf-8")

    print("Starting agent subprocess")
    a_proc.proc = await asyncio.create_subprocess_exec(
        sys.executable,
        (dir_p / "agent.py"),
        str(sock_agent_fd),
        pass_fds=[sock_agent_fd],
        stdin=asyncio.subprocess.DEVNULL,
        stdout=a_proc.log_io,
        stderr=a_proc.log_io,
        preexec_fn=lambda: os.nice(5),
    )
    async with a_proc.started:
        a_proc.started.notify_all()

    _ = Path(f"/proc/{a_proc.proc.pid}/oom_score_adj").write_text(
        "200\n", encoding="utf-8"
    )


async def bootstrap_headless_browser_on_startup() -> None:
    try:
        resp = await gql_query(
            auth=auth_token_sdk,
            query="""
                query AgentBootstrapPodInfo($podId: BigInt!) {
                    podInfos(filter: { id: { equalTo: $podId } }, first: 1) {
                        nodes {
                            plotNotebookId
                            accountId
                        }
                    }
                }
            """,
            variables={"podId": pod_id},
        )
        data = validate(resp, AgentBootstrapPodInfoResp).data
        nodes = data.podInfos.nodes
        if len(nodes) == 0:
            raise RuntimeError("No pod found for headless browser bootstrap")

        node = nodes[0]
        notebook_id = node.plotNotebookId
        workspace_id = node.accountId
        if notebook_id is None or workspace_id is None:
            raise RuntimeError("Missing notebook or workspace id for bootstrap")

        local_storage = {
            "plots.is_agent_controlled": "yes",
            "viewAccountId": str(workspace_id),
            "latch.authData": orjson.dumps({
                "status": "done",
                "auth0Data": {
                    "idToken": sdk_token.strip(),
                    "idTokenPayload": {
                        "sub": "agent-session",
                        "latch.bio/tos_ok": "true",
                    },
                },
            }).decode(),
        }

        print("[entrypoint] Starting headless browser during startup")
        await start_headless_browser(str(notebook_id), local_storage=local_storage)
    except Exception as e:
        print(f"[entrypoint] Failed to bootstrap headless browser: {e!s}")


async def start_headless_browser(
    notebook_id: str, local_storage: dict[str, str]
) -> None:
    global headless_browser, headless_browser_notebook_id, latest_local_storage

    local_storage_changed = local_storage and local_storage != latest_local_storage

    if not local_storage_changed and headless_browser is not None:
        return

    if local_storage_changed:
        latest_local_storage = local_storage
        print("[entrypoint] Local storage changed")

    if headless_browser is not None:
        print("[entrypoint] Stopping existing headless browser before starting new one")
        with contextlib.suppress(Exception):
            await headless_browser.stop()
        headless_browser = None

    headless_browser_notebook_id = notebook_id

    try:
        notebook_url = f"https://console.latch.bio/plots/{notebook_id}"

        headless_browser = HeadlessBrowser()
        await headless_browser.start(notebook_url, local_storage=local_storage)

        try:
            await asyncio.wait_for(action_handler_ready_ev.wait(), timeout=30)
            print(
                "[entrypoint] Headless browser action_handler connected successfully!"
            )
        except TimeoutError:
            print(
                "[entrypoint] Timed out waiting for action_handler websocket connection"
            )

    except Exception as e:
        print(f"[entrypoint] Error starting headless browser: {e}")
        traceback.print_exc()
        if headless_browser is not None:
            with contextlib.suppress(Exception):
                await headless_browser.screenshot("/var/log/headless_browser_error.png")
            with contextlib.suppress(Exception):
                await headless_browser.stop()
            headless_browser = None


async def restart_headless_browser() -> None:
    global headless_browser

    notebook_id = headless_browser_notebook_id
    local_storage = latest_local_storage

    if notebook_id is None or local_storage is None:
        print(
            "[entrypoint] Cannot restart headless browser: missing notebook_id or local_storage"
        )
        return

    print("[entrypoint] Restarting headless browser after disconnect...")

    if headless_browser is not None:
        with contextlib.suppress(Exception):
            await headless_browser.stop()
        headless_browser = None

    await start_headless_browser(notebook_id, local_storage)


async def poll_skills_branch() -> None:
    skills_dir = Path("/opt/latch/plots-faas/.claude/skills")
    skills_dir.mkdir(parents=True, exist_ok=True)

    skills_branch = "main"
    if sdk_token is not None:
        try:
            resp = await gql_query(
                auth=sdk_token,
                query="""
                    query GetNotebookMetadata($podId: BigInt!) {
                        podInfo(id: $podId) {
                            plotNotebook { metadata }
                        }
                    }
                """,
                variables={"podId": pod_id},
            )
            metadata_str = (
                resp
                .get("data", {})
                .get("podInfo", {})
                .get("plotNotebook", {})
                .get("metadata")
            )
            if metadata_str is not None:
                skills_branch = json.loads(metadata_str).get("skillsBranch", "main")
        except Exception as e:
            print(f"failed to fetch notebook metadata: {e}", file=sys.stderr)

    # todo: surface an error to the user if skills repo fails to pull or clone
    latch_skills_dest = skills_dir / "latch-skills"
    if latch_skills_dest.exists():
        ret = os.system(
            f"git -C {latch_skills_dest} fetch --depth 1 origin {skills_branch} && git -C {latch_skills_dest} checkout FETCH_HEAD"
        )
        if ret == 0:
            print(f"updated latch-skills to {skills_branch}")
        else:
            print(f"failed to update latch-skills to {skills_branch}", file=sys.stderr)
    else:
        ret = os.system(
            f"git clone --depth 1 --branch {skills_branch} https://github.com/latchbio/latch-skills.git {latch_skills_dest}"
        )
        if ret == 0:
            print(
                f"cloned public latch-skills ({skills_branch}) -> {latch_skills_dest}"
            )
        else:
            print("failed to clone public latch-skills", file=sys.stderr)

    if latch_skills_dest.exists():
        for link in skills_dir.iterdir():
            if link.is_symlink() and str(latch_skills_dest) in str(link.readlink()):
                link.unlink()

        for sub in sorted(latch_skills_dest.iterdir()):
            if sub.is_dir() and (sub / "SKILL.md").exists():
                link = skills_dir / sub.name
                if not link.exists():
                    link.symlink_to(sub)
                    print(f"linked latch skill: {sub.name} -> {link}")

    if sdk_token is not None:
        try:
            # todo(rteqs): combine into one query
            resp = await gql_query(
                auth=sdk_token,
                query="""
                    query AgentSkillRepos {
                        accountInfoCurrent {
                            id
                        }
                    }
                """,
                variables={},
            )
            account_id = resp.get("data", {}).get("accountInfoCurrent", {}).get("id")

            if account_id is None:
                raise RuntimeError("could not resolve account")

            resp = await gql_query(
                auth=sdk_token,
                query="""
                    query AgentSkillReposForAccount($accountId: BigInt!) {
                        agentSkillReposForAccount(argAccountId: $accountId)
                    }
                """,
                variables={"accountId": account_id},
            )
            skill_data = resp.get("data", {}).get("agentSkillReposForAccount")

            if skill_data is not None:
                pat = skill_data.get("pat")
                repos = skill_data.get("repos", [])

                if pat and repos:
                    username = pat["username"]
                    token = pat["token"]
                    skills_dir.mkdir(parents=True, exist_ok=True)

                    seen_dirs: set[str] = set()
                    for repo in repos:
                        repo_url: str = repo["repoUrl"]
                        authed_url = repo_url.replace(
                            "https://", f"https://{username}:{token}@"
                        )

                        parts = repo_url.rstrip("/").removesuffix(".git").split("/")
                        repo_name = parts[-1] if parts else repo["displayName"]

                        if repo_name in seen_dirs:
                            print(
                                f"skill repo conflict: {repo_name} already cloned, skipping {repo_url}",
                                file=sys.stderr,
                            )
                            continue
                        seen_dirs.add(repo_name)

                        dest = skills_dir / repo_name
                        if not dest.exists():
                            ret = os.system(f"git clone --depth 1 {authed_url} {dest}")
                            if ret == 0:
                                print(f"cloned skill repo: {repo_url} -> {dest}")
                            else:
                                print(
                                    f"failed to clone skill repo: {repo_url}",
                                    file=sys.stderr,
                                )
                                continue

                        if not (dest / "SKILL.md").exists():
                            for sub in sorted(dest.iterdir()):
                                if sub.is_dir() and (sub / "SKILL.md").exists():
                                    link = skills_dir / sub.name
                                    if sub.name in seen_dirs:
                                        print(
                                            f"skill repo conflict: {sub.name} already exists, skipping {sub}",
                                            file=sys.stderr,
                                        )
                                        continue
                                    seen_dirs.add(sub.name)
                                    if not link.exists():
                                        link.symlink_to(sub)
                                        print(
                                            f"linked monorepo skill: {sub.name} -> {link}"
                                        )
        except Exception as e:
            print(f"failed to fetch skill repos: {e}", file=sys.stderr)


async def start_poll_skills_branch() -> None:
    while True:
        await poll_skills_branch()
        await asyncio.sleep(60)


async def startup() -> None:
    await start_kernel_proc()

    async def kernel_died() -> None:
        for k, v in cell_status.items():
            if v != "running":
                continue

            cell_status[k] = "ran"

    async_tasks.append(asyncio.create_task(start_poll_skills_branch()))
    async_tasks.append(asyncio.create_task(k_proc.watchdog(cleanup=kernel_died)))
    async_tasks.append(asyncio.create_task(a_proc.watchdog()))
    async_tasks.append(asyncio.create_task(bootstrap_headless_browser_on_startup()))


async def shutdown() -> None:
    global shutting_down
    shutting_down = True

    await k_proc.stop()
    await a_proc.stop()

    for task in async_tasks:
        _ = task.cancel()

    _ = await asyncio.gather(*async_tasks)

    global headless_browser
    if headless_browser is not None:
        with contextlib.suppress(Exception):
            await headless_browser.stop()
        headless_browser = None

    with contextlib.suppress(Exception):
        sess = get_global_http_sess()
        await sess.close()
