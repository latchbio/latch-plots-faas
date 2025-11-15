import asyncio
import contextlib
import json
import os
import socket
import subprocess
import sys
import time
import traceback
import uuid
from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from agent_utils.auto_install import anthropic
from anthropic.types import MessageParam, ToolParam
from anthropic.types.beta.beta_message import BetaMessage
from anthropic.types.message import Message
from lplots import _inject
from socketio_thread import SocketIoThread
from utils import auth_token_sdk, gql_query, nucleus_url, pod_id

sys.stdout.reconfigure(line_buffering=True)

sandbox_root = os.environ.get("LATCH_SANDBOX_ROOT")
if sandbox_root:
    import pathlib
    original_path_new = pathlib.Path.__new__

    def patched_path_new(cls, *args, **kwargs):
        if args and args[0] == "/root/.latch":
            return original_path_new(cls, sandbox_root, *args[1:], **kwargs)
        return original_path_new(cls, *args, **kwargs)

    pathlib.Path.__new__ = patched_path_new


context_root = Path(__file__).parent / "agent_config" / "context"


class Mode(Enum):
    planning = "planning"
    executing = "executing"
    debugging = "debugging"


@dataclass
class AgentHarness:
    conn: SocketIoThread
    initialized: bool = False
    client: anthropic.AsyncAnthropic | None = None
    mode: Mode = Mode.planning
    system_prompt: str | None = None
    pending_operations: dict[str, asyncio.Future] = field(default_factory=dict)
    executing_cells: set[str] = field(default_factory=set)
    tools: list[ToolParam] = field(default_factory=list)
    tool_map: dict[str, Callable] = field(default_factory=dict)
    operation_counter: int = 0
    current_request_id: str | None = None
    should_auto_continue: bool = False
    pending_auto_continue: bool = False
    pending_tool_calls: set[str] = field(default_factory=set)

    pending_messages: asyncio.Queue = field(default_factory=asyncio.Queue)
    conversation_task: asyncio.Task | None = None
    conversation_running: bool = False
    agent_session_id: int | None = None
    latest_notebook_context: dict = field(default_factory=dict)
    current_status: str | None = None
    expected_widgets: dict[str, object | None] = field(default_factory=dict)

    mode_config: dict[Mode, tuple[str, int | None]] = field(default_factory=lambda: {
        Mode.planning: ("claude-sonnet-4-5-20250929", 4096),
        Mode.executing: ("claude-sonnet-4-5-20250929", 1024),
        Mode.debugging: ("claude-sonnet-4-5-20250929", 2048),
    })

    async def send(self, msg: dict[str, object]) -> None:
        msg_type = msg.get("type", "unknown")
        print(f"[agent] Sending message: {msg_type}")
        await self.conn.send(msg)

    async def _notify_history_updated(self, *, request_id: str | None = None) -> None:
        await self.send({
            "type": "agent_history_updated",
            "session_id": str(self.agent_session_id) if self.agent_session_id is not None else None,
            **({"request_id": request_id} if request_id else {}),
        })

    async def _fetch_history_from_db(self) -> list[dict]:
        assert self.agent_session_id is not None
        resp = await gql_query(
            auth=auth_token_sdk,
            query="""
                query AgentHistory($sessionId: BigInt!) {
                    agentHistories(condition: {sessionId: $sessionId, removed: false}, orderBy: ID_ASC) {
                        nodes { id payload requestId }
                    }
                }
            """,
            variables={"sessionId": str(self.agent_session_id)},
        )
        nodes = resp.get("data", {}).get("agentHistories", {}).get("nodes", [])
        return [{"payload": n.get("payload"), "request_id": n.get("requestId")} for n in nodes]

    def _truncate_old_messages(self, messages: list[MessageParam], keep_recent_turns: int = 3) -> list[MessageParam]:
        turn_boundaries = []
        for i, msg in enumerate(messages):
            if msg.get("role") == "user":
                turn_boundaries.append(i)

        if len(turn_boundaries) <= keep_recent_turns:
            return messages

        keep_from_index = turn_boundaries[-(keep_recent_turns)]

        truncated = []
        for i, msg in enumerate(messages):
            if i == 0 or i >= keep_from_index:
                truncated.append(msg)
            else:
                truncated.append(self._truncate_message_content(msg))

        return truncated

    @staticmethod
    def _truncate_message_content(msg: MessageParam) -> MessageParam:
        content = msg.get("content")

        if isinstance(content, str):
            if len(content) > 500:
                return {"role": msg["role"], "content": content[:500] + "...[truncated]"}
            return msg

        if isinstance(content, list):
            truncated_blocks = []
            for block in content:
                if not isinstance(block, dict):
                    truncated_blocks.append(block)
                    continue

                block_type = block.get("type")

                if block_type in {"thinking", "redacted_thinking"}:
                    continue

                if block_type == "tool_result":
                    result_str = block.get("content", "{}")
                    result = json.loads(result_str)
                    truncated_result = {
                        "tool_name": result.get("tool_name"),
                        "success": result.get("success"),
                        "summary": result.get("summary"),
                    }
                    truncated_blocks.append({
                        "type": "tool_result",
                        "tool_use_id": block.get("tool_use_id"),
                        "content": json.dumps(truncated_result)
                    })
                    continue

                if block_type == "tool_use":
                    truncated_block = block.copy()
                    if "input" in truncated_block:
                        inp = truncated_block["input"]
                        if isinstance(inp, dict) and "code" in inp and isinstance(inp["code"], str) and len(inp["code"]) > 1000:
                            truncated_block["input"] = {**inp, "code": inp["code"][:200] + "...[truncated]"}

                    truncated_blocks.append(truncated_block)
                    continue

                if block_type == "text":
                    text = block.get("text", "")
                    if len(text) > 1000:
                        truncated_blocks.append({"type": "text", "text": text[:1000] + "...[truncated]"})
                    else:
                        truncated_blocks.append(block)
                    continue

                truncated_blocks.append(block)

            return {"role": msg["role"], "content": truncated_blocks}

        return msg

    async def _build_messages_from_db(self) -> list[MessageParam]:
        history = await self._fetch_history_from_db()
        anthropic_messages: list[MessageParam] = []

        for item in history:
            payload = item.get("payload")

            t = payload.get("type") if isinstance(payload, dict) else None
            if t == "anthropic_message":
                role = payload.get("role")
                content = payload.get("content")

                if (role == "user" and isinstance(content, dict) and content.get("type") == "cell_result"):
                    exception = content.get("exception")
                    logs = content.get("logs")
                    message = content.get("message", "Cell execution completed")

                    content = message
                    if exception:
                        content = f"{message}\n\nException: {exception}"
                    if logs:
                        content = f"{content}\n\nLogs:\n{logs}"

                if isinstance(content, list):
                    cleaned_content = []
                    for block in content:
                        if isinstance(block, dict) and block.get("type") == "tool_result":
                            block = block.copy()
                            result = json.loads(block.get("content", "{}"))
                            if "original_code" in result:
                                result.pop("original_code")
                                block["content"] = json.dumps(result)
                        cleaned_content.append(block)
                    content = cleaned_content

                if role in {"user", "assistant"} and (isinstance(content, (str, list))):
                    anthropic_messages.append({"role": role, "content": content})

        print(f"[agent] Built {len(anthropic_messages)} messages from DB")

        return anthropic_messages

    async def _extract_last_request_id(self) -> str | None:
        history = await self._fetch_history_from_db()

        for item in reversed(history):
            if not isinstance(item, dict):
                continue
            request_id = item.get("request_id")
            if request_id is not None:
                return request_id

        return None

    async def _insert_history(
        self,
        *,
        event_type: str = "anthropic_message",
        role: str = "user",
        payload: dict,
        request_id: str | None = None,
        tx_id: str | None = None,
    ) -> None:
        assert self.agent_session_id is not None

        payload = {
            "type": event_type,
            "role": role,
            "timestamp": int(time.time() * 1000),
            **payload,
        }

        variables = {
            "sessionId": str(self.agent_session_id),
            "eventType": event_type,
            "payload": payload,
            "requestId": request_id,
            "txId": tx_id,
        }

        await gql_query(
            auth=auth_token_sdk,
            query="""
                mutation CreateAgentHistory($sessionId: BigInt!, $eventType: String!, $payload: JSON!, $requestId: String, $txId: String) {
                    createAgentHistory(input: {agentHistory: {sessionId: $sessionId, eventType: $eventType, payload: $payload, requestId: $requestId, txId: $txId}}) {
                        clientMutationId
                    }
                }
            """,
            variables=variables,
        )

        await self._notify_history_updated(request_id=request_id)

    async def _mark_all_history_removed(self) -> None:
        assert self.agent_session_id is not None

        await gql_query(
            auth=auth_token_sdk,
            query="""
                mutation ClearAgentHistory($sessionId: BigInt!) {
                    clearAgentHistory(input: {argSessionId: $sessionId}) {
                        clientMutationId
                    }
                }
            """,
            variables={"sessionId": str(self.agent_session_id)},
        )
        await self._notify_history_updated()

    def _start_conversation_loop(self) -> None:
        self.conversation_task = asyncio.create_task(self.run_agent_loop())

        def _task_done_callback(task: asyncio.Task) -> None:
            try:
                task.result()
            except Exception as e:
                print(f"[agent] conversation_task raised exception: {e}")
                traceback.print_exc()

        self.conversation_task.add_done_callback(_task_done_callback)

    async def _clear_running_state(self) -> None:
        if len(self.pending_operations) > 0:
            print(f"[agent] Cancelling {len(self.pending_operations)} pending operations")
            for tx_id, future in self.pending_operations.items():
                if not future.done():
                    future.cancel()
                    print(f"[agent]   Cancelled: {tx_id}")
            self.pending_operations.clear()

        if self.conversation_task and not self.conversation_task.done():
            print("[agent] Cancelling conversation task")
            self.conversation_running = False
            self.conversation_task.cancel()
            with contextlib.suppress(asyncio.CancelledError):
                await self.conversation_task
            self.conversation_task = None

        while not self.pending_messages.empty():
            try:
                self.pending_messages.get_nowait()
            except asyncio.QueueEmpty:
                break

        self.executing_cells.clear()

        print("[agent] Cleared running state")

    async def atomic_operation(self, action: str, params: dict | None = None) -> dict:
        if params is None:
            params: dict = {}

        self.operation_counter += 1

        if self.mode == Mode.planning and action in {"create_cell", "edit_cell", "run_cell", "delete_cell"}:
            self.set_mode(Mode.executing)

        tx_id = f"tx_{uuid.uuid4().hex[:12]}"
        loop = asyncio.get_running_loop()
        response_future = loop.create_future()
        self.pending_operations[tx_id] = response_future

        start_time = time.time()
        try:
            print(f"[agent] -> {action}")
            await self.send({"type": "agent_action", "action": action, "params": params, "tx_id": tx_id})
        except Exception as e:
            self.pending_operations.pop(tx_id, None)
            return {"status": "error", "error": f"Send failed: {e!s}"}

        try:
            return await asyncio.wait_for(response_future, timeout=10.0)
        except asyncio.CancelledError:
            print(f"[agent] Operation cancelled (session reinitialized): action={action}, tx_id={tx_id}")
            return {"status": "error", "error": f"OPERATION FAILED: '{action}' was interrupted because the session was reinitialized. This operation did NOT complete. You must retry or inform the user."}
        except TimeoutError:
            return {"status": "error", "error": f"OPERATION FAILED: '{action}' timed out after 10 seconds. This operation did NOT complete.", "tx_id": tx_id}
        finally:
            duration = time.time() - start_time
            print(f"[agent] {action} took {duration:.3f}s")
            self.pending_operations.pop(tx_id, None)

    async def handle_action_response(self, msg: dict[str, object]) -> None:
        tx_id = msg.get("tx_id")
        fut = self.pending_operations.get(tx_id)
        if fut and not fut.done():
            fut.set_result(msg)

    def set_mode(self, mode: Mode) -> None:
        if mode == self.mode:
            return

        self.mode = mode
        print(f"[agent] Mode changed to {mode.value}")

    async def _wait_for_message(self) -> None:
        print("[agent] _wait_for_message: waiting for message...")
        msg = await self.pending_messages.get()

        msg_type = msg.get("type")

        if msg_type == "resume":
            print("[agent] Resuming turn after tool results")
            return

        if msg_type == "user_query":
            self.current_request_id = msg.get("request_id")

            payload = {
                "content": msg["content"],
            }

            if msg.get("display_query") is not None:
                payload["display_query"] = msg["display_query"]
            if msg.get("display_nodes") is not None:
                payload["display_nodes"] = msg["display_nodes"]
            if msg.get("hidden") is not None:
                payload["hidden"] = msg["hidden"]

            await self._insert_history(
                payload=payload,
                request_id=self.current_request_id,
            )

        elif msg_type == "cell_result":
            if self.current_request_id is None:
                print(f"[agent] Ignoring cell_result for {msg.get('cell_id')} - no active request")
                return

            cell_id = msg["cell_id"]
            success = msg.get("success", True)
            cell_name = msg.get("display_name", None)

            if success:
                result_message = f"✓ Cell {cell_name} ({cell_id}) executed successfully"
                result_content = {
                    "type": "cell_result",
                    "message": result_message,
                    "cell_id": cell_id,
                    "cell_name": cell_name,
                    "success": success,
                    "logs": msg.get("logs"),
                }
                print(f"[agent] Cell {cell_id} succeeded")
            else:
                exception = msg.get("exception", "Unknown error")
                result_message = f"✗ Cell {cell_name} ({cell_id}) execution failed"
                result_content = {
                    "type": "cell_result",
                    "message": result_message,
                    "cell_id": cell_id,
                    "cell_name": cell_name,
                    "success": False,
                    "exception": exception,
                    "logs": msg.get("logs"),
                }
                print(f"[agent] Cell {cell_id} failed")

            await self._insert_history(
                payload={
                    "content": result_content,
                },
            )

            if self.pending_auto_continue and not self.executing_cells:
                print(f"[agent] All cells complete, resuming auto-continue (request_id={self.current_request_id})")
                self.pending_auto_continue = False
                await self.pending_messages.put({
                    "type": "user_query",
                    "content": "Continue with the next step.",
                    "request_id": self.current_request_id,
                    "hidden": True,
                })

        elif msg_type == "set_widget_value":
            data = msg.get("data", {})
            widget_info = ", ".join(f"{k}={v}" for k, v in data.items())
            content = f"User provided input via widget(s): {widget_info}"

            await self._insert_history(
                payload={
                    "content": content,
                    "hidden": True,
                },
            )

    async def _complete_turn(self) -> None:
        if self.current_request_id is None:
            return

        should_continue = self.should_auto_continue
        self.should_auto_continue = False
        self.pending_auto_continue = False

        if self.mode == Mode.executing:
            self.set_mode(Mode.planning)

        await self._notify_history_updated(request_id=self.current_request_id)

        if should_continue:
            print("[agent] Auto-continuing as requested by model")
            await self.pending_messages.put({
                "type": "user_query",
                "content": "Continue with the next step.",
                "request_id": self.current_request_id,
                "hidden": True,
            })

    def init_tools(self) -> None:
        self.tools = []
        self.tool_map = {}

        async def create_cell(args: dict) -> dict:
            position = args["position"]
            code = args["code"]
            title = args["title"]
            action_summary = args["action_summary"]

            if position < 0:
                return {
                    "tool_name": "create_cell",
                    "summary": "Error: Position must be non-negative",
                    "success": False,
                }

            print(f'[tool] create_cell pos={position} title="{title}"')

            params = {
                "position": position,
                "cell_type": "code",
                "source": code,
                "title": title,
                "auto_run": True,
            }

            result = await self.atomic_operation("create_cell", params)
            if result.get("status") == "success":
                cell_id = result.get("cell_id", "unknown")
                msg = f"Created cell at position {position} (ID: {cell_id}, Title: {title})"
                print(f"[tool] create_cell -> {msg}")
                return {
                    "tool_name": "create_cell",
                    "summary": msg,
                    "code": code,
                    "cell_id": cell_id,
                    "cell_name": title,
                    "position": position,
                    "message": action_summary,
                    "success": True,
                }
            return {
                "tool_name": "create_cell",
                "summary": f"Failed to create cell: {result.get('error', 'Unknown error')}",
                "success": False,
            }

        async def create_markdown_cell(args: dict) -> dict:
            position = args["position"]
            code = args["code"]
            title = args["title"]
            action_summary = args["action_summary"]

            if position < 0:
                return {
                    "summary": "Error: Position must be non-negative",
                    "success": False,
                }

            print(f"[tool] create_markdown_cell pos={position}")

            params = {
                "position": position,
                "cell_type": "markdown",
                "source": code,
            }

            result = await self.atomic_operation("create_markdown_cell", params)
            if result.get("status") == "success":
                cell_id = result.get("cell_id", "unknown")
                msg = f"Created markdown cell at position {position} (ID: {cell_id})"

                print(f"[tool] create_markdown_cell -> {msg}")
                return {
                    "tool_name": "create_markdown_cell",
                    "summary": msg,
                    "code": code,
                    "cell_id": cell_id,
                    "cell_name": title,
                    "position": position,
                    "message": action_summary,
                    "success": True,
                }
            return {
                "tool_name": "create_markdown_cell",
                "message": f"Failed to create cell: {result.get('error', 'Unknown error')}",
                "success": False,
            }

        async def edit_cell(args: dict) -> dict:
            cell_id = args["cell_id"]
            new_code = args["new_code"]
            title = args["title"]
            action_summary = args["action_summary"]

            print(f"[tool] edit_cell id={cell_id}")

            original_code = None

            # ensure self.latest_notebook_context is up to date
            await self.tool_map.get("refresh_cells_context")({})
            cells = self.latest_notebook_context.get("cells", [])

            for cell in cells:
                if cell.get("cell_id") == cell_id:
                    original_code = cell.get("source", "")
                    break

            params = {
                "cell_id": cell_id,
                "source": new_code,
                "auto_run": True
            }

            result = await self.atomic_operation("edit_cell", params)
            if result.get("status") == "success":
                msg = f"Cell {cell_id} edited successfully"
                print(f"[tool] edit_cell -> {msg}")

                return {
                    "tool_name": "edit_cell",
                    "summary": msg,
                    "code": new_code,
                    "original_code": original_code,
                    "cell_id": cell_id,
                    "cell_name": title,
                    "message": action_summary,
                    "success": True,
                }
            return {
                "tool_name": "edit_cell",
                "summary": f"Failed to edit cell: {result.get('error', 'Unknown error')}",
                "success": False,
            }

        async def delete_cell(args: dict) -> dict:
            cell_id = args["cell_id"]
            title = args["title"]
            action_summary = args["action_summary"]

            print(f"[tool] delete_cell id={cell_id}")
            params = {"cell_id": cell_id}

            result = await self.atomic_operation("delete_cell", params)
            if result.get("status") == "success":
                remaining = result.get("remaining_cells", [])
                cell_count = result.get("cell_count", 0)

                if remaining:
                    cell_list = ", ".join([f"{c['index']}: {c['cell_type']}" for c in remaining[:5]])
                    if len(remaining) > 5:
                        cell_list += f", ... ({len(remaining) - 5} more)"
                    msg = f"Cell {cell_id} deleted. {cell_count} cells remain: [{cell_list}]"
                else:
                    msg = f"Cell {cell_id} deleted. No cells remain in notebook."
                print(f"[tool] delete_cell -> {msg}")
                return {
                    "tool_name": "delete_cell",
                    "summary": msg,
                    "cell_id": cell_id,
                    "cell_name": title,
                    "message": action_summary,
                    "success": True,
                }
            return {
                "tool_name": "delete_cell",
                "summary": f"Failed to delete cell: {result.get('error', 'Unknown error')}",
                "success": False,
            }

        async def run_cell(args: dict) -> dict:
            cell_id = args["cell_id"]
            title = args["title"]
            action_summary = args["action_summary"]
            params = {"cell_id": cell_id}

            await self.send({
                "type": "agent_action",
                "action": "run_cell",
                "params": params,
            })
            self.executing_cells.add(cell_id)

            return {
                "tool_name": "run_cell",
                "summary": f"Cell {cell_id} execution started",
                "cell_id": cell_id,
                "cell_name": title,
                "message": action_summary,
                "success": True,
            }

        async def stop_cell(args: dict) -> dict:
            cell_id = args["cell_id"]
            title = args["title"]
            action_summary = args["action_summary"]
            params = {"cell_id": cell_id}

            result = await self.atomic_operation("stop_cell", params)
            if result.get("status") == "success":
                self.executing_cells.discard(cell_id)
                return {
                    "tool_name": "stop_cell",
                    "summary": f"Stopped cell {cell_id}",
                    "cell_id": cell_id,
                    "cell_name": title,
                    "message": action_summary,
                }
            return {
                "tool_name": "stop_cell",
                "summary": f"Failed to stop cell {cell_id}: {result.get('error', 'Unknown error')}",
                "success": False,
            }

        async def delete_all_cells(args: dict) -> dict:
            context_result = await self.atomic_operation("get_context", {})
            if context_result.get("status") != "success":
                error_msg = context_result.get("error", "Unknown error")
                return {
                    "tool_name": "delete_all_cells",
                    "summary": f"Failed to delete cells: {error_msg}",
                    "success": False,
                }

            cells = context_result.get("context", {}).get("cells", [])
            deleted_count = 0

            for cell in reversed(cells):
                cell_id = cell.get("cell_id")
                if cell_id:
                    result = await self.atomic_operation("delete_cell", {"cell_id": cell_id})
                    if result.get("status") == "success":
                        deleted_count += 1

            return {
                "tool_name": "delete_all_cells",
                "success": True,
                "summary": f"Deleted {deleted_count} cells from the notebook",
                "deleted_count": deleted_count,
            }

        async def rename_notebook(args: dict) -> dict:
            name = args["name"]
            print(f"[tool] rename_notebook name={name}")
            params = {"name": name}

            result = await self.atomic_operation("rename_notebook", params)
            if result.get("status") == "success":
                return {
                    "tool_name": "rename_notebook",
                    "success": True,
                    "summary": f"Notebook renamed to '{name}'",
                    "name": name,
                }

            return {
                "tool_name": "rename_notebook",
                "summary": f"Failed to rename notebook: {result.get('error', 'Unknown error')}",
                "success": False,
            }

        async def create_tab(args: dict) -> dict:
            position = args["position"]
            display_name = args["display_name"]

            if position < 0:
                return {
                    "tool_name": "create_tab",
                    "summary": "Error: Position must be non-negative",
                    "success": False,
                }

            print(f'[tool] create_tab pos={position} name="{display_name}"')

            params = {
                "position": position,
                "display_name": display_name,
            }

            result = await self.atomic_operation("create_tab", params)
            if result.get("status") == "success":
                tab_id = result.get("tab_id", "unknown")
                msg = f"Created tab at position {position} (ID: {tab_id}, Name: {display_name})"
                print(f"[tool] create_tab -> {msg}")
                return {
                    "tool_name": "create_tab",
                    "summary": msg,
                    "tab_id": tab_id,
                    "display_name": display_name,
                    "position": position,
                    "success": True,
                }
            return {
                "tool_name": "create_tab",
                "summary": f"Failed to create tab: {result.get('error', 'Unknown error')}",
                "success": False,
            }

        async def rename_tab(args: dict) -> dict:
            tab_id = args["tab_id"]
            new_name = args["new_name"]

            print(f'[tool] rename_tab tab_id={tab_id} new_name="{new_name}"')

            params = {
                "tab_id": tab_id,
                "new_name": new_name,
            }

            result = await self.atomic_operation("rename_tab", params)
            if result.get("status") == "success":
                target = "Tab 1" if tab_id == "DEFAULT" else f"tab {tab_id}"
                msg = f"Renamed {target} to '{new_name}'"
                print(f"[tool] rename_tab -> {msg}")
                return {
                    "tool_name": "rename_tab",
                    "summary": msg,
                    "tab_id": tab_id,
                    "new_name": new_name,
                    "success": True,
                }
            return {
                "tool_name": "rename_tab",
                "summary": f"Failed to rename tab: {result.get('error', 'Unknown error')}",
                "success": False,
            }

        async def set_widget(args: dict) -> dict:
            key = args.get("key")
            action_summary = args.get("action_summary")
            label = args.get("label")
            if not key:
                return {
                    "tool_name": "set_widget",
                    "summary": "Failed to set widget: Widget key is required",
                    "success": False,
                }

            value = args.get("value")
            if value is None:
                return {
                    "tool_name": "set_widget",
                    "summary": "Failed to set widget: Widget value is required",
                    "success": False,
                }

            print(f"[tool] set_widget key={key} value={value!r}")

            params = {"key": key, "value": json.dumps(value)}

            result = await self.atomic_operation("set_widget", params)

            if result.get("status") == "success":
                return {
                    "tool_name": "set_widget",
                    "success": True,
                    "summary": f"Updated widget value for: {key}",
                    "key": key,
                    "label": label,
                    "value": value,
                    "message": action_summary,
                }
            return {
                "tool_name": "set_widget",
                "summary": f"Failed to update widget value: {result.get('error', 'Unknown error')}",
                "success": False,
            }

        def submit_response(args: dict) -> dict:
            try:
                summary = args.get("summary")
                if summary is not None and not isinstance(summary, str):
                    summary = None

                questions = args.get("questions")
                if questions is not None and not isinstance(questions, str):
                    questions = None

                next_status = args.get("next_status")
                if not isinstance(next_status, str) or next_status not in {"executing", "fixing", "thinking", "awaiting_user_response", "awaiting_cell_execution", "awaiting_user_widget_input", "done"}:
                    print(f"[agent] Invalid next_status: {next_status}")
                    return {
                        "tool_name": "submit_response",
                        "message": "Please provide a valid next_status",
                        "success": False,
                    }

                should_continue = args.get("continue", False)

                plan_items = args.get("plan", [])
                plan_diff_items = args.get("plan_diff", [])
                expected_widgets = args.get("expected_widgets", [])

                print("[tool] submit_response called with:")
                print(f"  - next_status: {next_status}")
                print(f"  - plan: {len(plan_items)} items")
                print(f"  - plan_update_overview: {args.get('plan_update_overview')}")
                for item in plan_items:
                    print(f"    - [{item.get('status')}] {item.get('id')}: {item.get('description')}")
                print(f"  - plan_diff: {len(plan_diff_items)} items")
                for diff in plan_diff_items:
                    print(f"    - [{diff.get('action')}] {diff.get('id')}: {diff.get('description')}")
                print(f"  - summary: {summary}")
                print(f"  - questions: {questions}")
                print(f"  - continue: {should_continue}")
                print(f"  - expected_widgets: {expected_widgets}")

                if should_continue and self.executing_cells:
                    print(f"[tool] Deferring auto-continue - {len(self.executing_cells)} cells still executing: {self.executing_cells}")
                    self.should_auto_continue = False
                    self.pending_auto_continue = True
                else:
                    self.should_auto_continue = should_continue
                    self.pending_auto_continue = False

                self.current_status = next_status
                if next_status == "awaiting_user_widget_input":
                    self.expected_widgets = {str(k): None for k in expected_widgets}

                return {
                    "tool_name": "submit_response",
                    "summary": "Response submitted successfully",
                    "success": True,
                }
            except Exception as e:
                print(f"[tool] submit_response error: {e}")
                import traceback
                traceback.print_exc()
                return {
                    "tool_name": "submit_response",
                    "summary": f"Error submitting response: {e!s}",
                    "success": False,
                }

        async def h5_filter_by(args: dict) -> dict:
            widget_key = args.get("widget_key")
            filters = args.get("filters")

            if isinstance(filters, str):
                try:
                    filters = json.loads(filters)
                except json.JSONDecodeError:
                    return {
                        "tool_name": "h5_filter_by",
                        "success": False,
                        "summary": f"filters is invalid JSON: {filters!r}",
                    }

            print(f"[tool] h5_filter_by widget_key={widget_key} filters={filters}")

            params = {
                "widget_key": widget_key,
                "filters": filters
            }

            result = await self.atomic_operation("h5_filter_by", params)
            if result.get("status") == "success":
                return {
                    "tool_name": "h5_filter_by",
                    "success": True,
                    "label": args.get("label"),
                    "summary": f"Applied filters to h5 widget: {filters}",
                    "widget_key": widget_key,
                    "filters": filters,
                }

            return {
                "tool_name": "h5_filter_by",
                "success": False,
                "summary": f"Failed to apply filters to h5 widget: {result.get('error', 'Unknown error')}",
            }

        async def h5_color_by(args: dict) -> dict:
            widget_key = args.get("widget_key")
            color_by = args.get("color_by")

            if isinstance(color_by, str):
                try:
                    color_by = json.loads(color_by)
                except json.JSONDecodeError:
                    return {
                        "tool_name": "h5_color_by",
                        "success": False,
                        "summary": f"color_by is invalid JSON: {color_by!r}",
                    }

            print(f"[tool] h5_color_by widget_key={widget_key} color_by={color_by}")

            params = {
                "widget_key": widget_key,
                "color_by": color_by
            }

            result = await self.atomic_operation("h5_color_by", params)
            if result.get("status") == "success":
                return {
                    "tool_name": "h5_color_by",
                    "success": True,
                    "label": args.get("label"),
                    "summary": f"Set h5 widget coloring: {color_by}",
                    "widget_key": widget_key,
                    "color_by": color_by,
                }

            return {
                "tool_name": "h5_color_by",
                "success": False,
                "summary": f"Failed to set h5 widget coloring: {result.get('error', 'Unknown error')}",
            }

        async def h5_set_selected_obsm_key(args: dict) -> dict:
            widget_key = args.get("widget_key")
            obsm_key = args.get("obsm_key")

            print(f"[tool] h5_set_selected_obsm_key widget_key={widget_key} obsm_key={obsm_key}")

            params = {
                "widget_key": widget_key,
                "obsm_key": obsm_key
            }

            result = await self.atomic_operation("h5_set_selected_obsm_key", params)
            if result.get("status") == "success":
                return {
                    "tool_name": "h5_set_selected_obsm_key",
                    "success": True,
                    "label": args.get("label"),
                    "summary": f"Set h5 widget to use obsm key {obsm_key}",
                    "widget_key": widget_key,
                    "obsm_key": obsm_key,
                }

            return {
                "tool_name": "h5_set_selected_obsm_key",
                "success": False,
                "summary": f"Failed to set h5 widget obsm key: {result.get('error', 'Unknown error')}",
            }

        async def h5_set_background_image(args: dict) -> dict:
            widget_key = args.get("widget_key")
            node_id = args.get("node_id")

            print(f"[tool] h5_set_background_image widget_key={widget_key} node_id={node_id}")

            params = {
                "widget_key": widget_key,
                "node_id": node_id,
            }

            result = await self.atomic_operation("h5_set_background_image", params)
            if result.get("status") == "success":
                return {
                    "tool_name": "h5_set_background_image",
                    "success": True,
                    "label": args.get("label"),
                    "summary": f"Set background image for h5 widget using {node_id}",
                    "widget_key": widget_key,
                    "node_id": node_id,
                }
            return {
                "tool_name": "h5_set_background_image",
                "success": False,
                "summary": f"Failed to set background image: {result.get('error', 'Unknown error')}",
            }

        async def h5_open_image_aligner(args: dict) -> dict:
            widget_key = args.get("widget_key")
            background_image_id = args.get("background_image_id")

            print(f"[tool] h5_open_image_aligner widget_key={widget_key} background_image_id={background_image_id}")

            params = {
                "widget_key": widget_key,
                "background_image_id": background_image_id,
            }

            result = await self.atomic_operation("h5_open_image_aligner", params)
            if result.get("status") == "success":
                return {
                    "tool_name": "h5_open_image_aligner",
                    "success": True,
                    "label": args.get("label"),
                    "summary": f"Opened image aligner for background image {background_image_id}",
                    "widget_key": widget_key,
                    "background_image_id": background_image_id,
                }

            return {
                "tool_name": "h5_open_image_aligner",
                "success": False,
                "summary": f"Failed to open image aligner: {result.get('error', 'Unknown error')}",
            }

        async def h5_autoscale(args: dict) -> dict:
            widget_key = args.get("widget_key")

            print(f"[tool] h5_autoscale widget_key={widget_key}")

            params = {
                "widget_key": widget_key,
            }

            result = await self.atomic_operation("h5_autoscale", params)
            if result.get("status") == "success":
                return {
                    "tool_name": "h5_autoscale",
                    "success": True,
                    "label": args.get("label"),
                    "summary": f"Autoscaled h5 widget {widget_key} to data bounds",
                    "widget_key": widget_key,
                }

            return {
                "tool_name": "h5_autoscale",
                "success": False,
                "summary": f"Failed to autoscale h5 widget: {result.get('error', 'Unknown error')}",
            }

        async def h5_zoom(args: dict) -> dict:
            widget_key = args.get("widget_key")
            direction = args.get("direction")
            percentage = args.get("percentage")

            print(f"[tool] h5_zoom widget_key={widget_key} direction={direction} percentage={percentage}")

            params = {
                "widget_key": widget_key,
                "direction": direction,
            }

            if percentage is not None:
                params["percentage"] = percentage

            result = await self.atomic_operation("h5_zoom", params)
            if result.get("status") == "success":
                zoom_desc = f"zoom {direction}"
                if percentage is not None:
                    zoom_desc += f" by {percentage}%"
                return {
                    "tool_name": "h5_zoom",
                    "success": True,
                    "label": args.get("label"),
                    "summary": f"Applied {zoom_desc} to h5 widget {widget_key}",
                    "widget_key": widget_key,
                    "direction": direction,
                    "percentage": percentage,
                }

            return {
                "tool_name": "h5_zoom",
                "success": False,
                "summary": f"Failed to zoom h5 widget: {result.get('error', 'Unknown error')}",
            }

        async def h5_set_background_image_visibility(args: dict) -> dict:
            widget_key = args.get("widget_key")
            background_image_id = args.get("background_image_id")
            hidden = args.get("hidden")

            print(f"[tool] h5_set_background_image_visibility widget_key={widget_key} background_image_id={background_image_id} hidden={hidden}")

            params = {
                "widget_key": widget_key,
                "background_image_id": background_image_id,
                "hidden": hidden
            }

            result = await self.atomic_operation("h5_set_background_image_visibility", params)
            if result.get("status") == "success":
                visibility_action = "hidden" if hidden else "shown"
                return {
                    "tool_name": "h5_set_background_image_visibility",
                    "success": True,
                    "label": args.get("label"),
                    "summary": f"Background image {background_image_id} {visibility_action}",
                    "widget_key": widget_key,
                    "background_image_id": background_image_id,
                    "hidden": hidden,
                }

            return {
                "tool_name": "h5_set_background_image_visibility",
                "success": False,
                "summary": f"Failed to set background image visibility: {result.get('error', 'Unknown error')}",
            }

        async def h5_add_selected_cells_to_categorical_obs(args: dict) -> dict:
            widget_key = args.get("widget_key")
            obs_key = args.get("obs_key")
            category = args.get("category")

            print(f"[tool] h5_add_selected_cells_to_categorical_obs widget_key={widget_key} obs_key={obs_key} category={category}")

            params = {
                "widget_key": widget_key,
                "obs_key": obs_key,
                "category": category
            }

            result = await self.atomic_operation("h5_add_selected_cells_to_categorical_obs", params)
            if result.get("status") == "success":
                return {
                    "tool_name": "h5_add_selected_cells_to_categorical_obs",
                    "success": True,
                    "label": args.get("label"),
                    "summary": f"Assigned selected cells to category '{category}' in observation key '{obs_key}'",
                    "widget_key": widget_key,
                    "obs_key": obs_key,
                    "category": category,
                }

            return {
                "tool_name": "h5_add_selected_cells_to_categorical_obs",
                "success": False,
                "summary": f"Failed to assign selected cells to category: {result.get('error', 'Unknown error')}",
            }

        async def h5_set_marker_opacity(args: dict) -> dict:
            widget_key = args.get("widget_key")
            opacity = args.get("opacity")

            print(f"[tool] h5_set_marker_opacity widget_key={widget_key} opacity={opacity}")

            params = {
                "widget_key": widget_key,
                "opacity": opacity
            }

            result = await self.atomic_operation("h5_set_marker_opacity", params)
            if result.get("status") == "success":
                return {
                    "tool_name": "h5_set_marker_opacity",
                    "success": True,
                    "label": args.get("label"),
                    "summary": f"Set marker opacity to {opacity}",
                    "widget_key": widget_key,
                    "opacity": opacity,
                }

            return {
                "tool_name": "h5_set_marker_opacity",
                "success": False,
                "summary": f"Failed to set marker opacity: {result.get('error', 'Unknown error')}",
            }

        async def h5_manage_obs(args: dict) -> dict:
            widget_key = args.get("widget_key")
            obs_key = args.get("obs_key")
            operation = args.get("operation")
            obs_type = args.get("obs_type", "category")

            print(f"[tool] h5_manage_obs widget_key={widget_key} obs_key={obs_key} operation={operation} obs_type={obs_type}")

            params = {
                "widget_key": widget_key,
                "obs_key": obs_key,
                "operation": operation
            }

            if operation == "add":
                params["obs_type"] = obs_type

            result = await self.atomic_operation("h5_manage_obs", params)
            if result.get("status") == "success":
                if operation == "add":
                    return {
                        "tool_name": "h5_manage_obs",
                        "success": True,
                        "label": args.get("label"),
                        "summary": f"Created observation column '{obs_key}' with type '{obs_type}'",
                        "widget_key": widget_key,
                        "obs_key": obs_key,
                        "operation": operation,
                        "obs_type": obs_type,
                    }
                else:
                    return {
                        "tool_name": "h5_manage_obs",
                        "success": True,
                        "label": args.get("label"),
                        "summary": f"Deleted observation column '{obs_key}'",
                        "widget_key": widget_key,
                        "obs_key": obs_key,
                        "operation": operation,
                    }

            return {
                "tool_name": "h5_manage_obs",
                "success": False,
                "summary": f"Failed to {operation} observation column: {result.get('error', 'Unknown error')}",
            }

        async def smart_ui_spotlight(args: dict) -> dict:
            keyword = args.get("keyword")
            widget_key = args.get("widget_key")

            print(f"[tool] smart_ui_spotlight keyword={keyword}, widget_key={widget_key}")

            params = {
                "keyword": keyword
            }
            if widget_key is not None:
                params["widget_key"] = widget_key

            result = await self.atomic_operation("smart_ui_spotlight", params)
            if result.get("status") == "success":
                return {
                    "tool_name": "smart_ui_spotlight",
                    "success": True,
                    "summary": f"Highlighted UI element: {keyword}",
                    "keyword": keyword,
                }

            return {
                "tool_name": "smart_ui_spotlight",
                "success": False,
                "summary": f"Failed to highlight UI element: {result.get('error', 'Unknown error')}",
            }

        self.tools.append({
            "name": "create_cell",
            "description": "Create a new code cell at specified position. The cell will automatically run after creation.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "position": {"type": "integer", "description": "Position to insert the cell"},
                    "code": {"type": "string", "description": "Python code for the cell"},
                    "title": {"type": "string", "description": "Name for the cell"},
                    "action_summary": {"type": "string", "description": "Summary of the purpose of the cell."},
                },
                "required": ["position", "code", "title", "action_summary"],
            },
        })
        self.tool_map["create_cell"] = create_cell

        self.tools.append({
            "name": "create_markdown_cell",
            "description": "Create a new markdown cell at specified position.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "position": {"type": "integer", "description": "Position to insert the cell"},
                    "code": {"type": "string", "description": "Markdown content"},
                    "title": {"type": "string", "description": "Title of first header in the markdown cell"},
                    "action_summary": {"type": "string", "description": "Summary of the purpose of the cell."},
                },
                "required": ["position", "code", "title", "action_summary"],
            },
        })
        self.tool_map["create_markdown_cell"] = create_markdown_cell

        self.tools.append({
            "name": "edit_cell",
            "description": "Replace the contents of an existing cell. The cell will automatically run after editing.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "cell_id": {"type": "string", "description": "ID of the cell to edit"},
                    "new_code": {"type": "string", "description": "New code/content for the cell"},
                    "title": {"type": "string", "description": "Name of the cell to edit"},
                    "action_summary": {"type": "string", "description": "Summary of the purpose of the edit."},
                },
                "required": ["cell_id", "new_code", "title"],
            },
        })
        self.tool_map["edit_cell"] = edit_cell

        self.tools.append({
            "name": "delete_cell",
            "description": "Remove a cell from the notebook.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "cell_id": {"type": "string", "description": "ID of the cell to delete"},
                    "title": {"type": "string", "description": "Name of the cell to delete"},
                    "action_summary": {"type": "string", "description": "Summary of the purpose of the delete."},
                },
                "required": ["cell_id", "title", "action_summary"],
            },
        })
        self.tool_map["delete_cell"] = delete_cell

        self.tools.append({
            "name": "run_cell",
            "description": "Execute a specific cell.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "cell_id": {"type": "string", "description": "ID of the cell to run"},
                    "title": {"type": "string", "description": "Name of the cell to run"},
                    "action_summary": {"type": "string", "description": "Summary of the purpose of the run."},
                },
                "required": ["cell_id", "title", "action_summary"],
            },
        })
        self.tool_map["run_cell"] = run_cell

        self.tools.append({
            "name": "stop_cell",
            "description": "Stop execution of a specific cell.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "cell_id": {"type": "string", "description": "ID of the cell to stop"},
                    "cell_name": {"type": "string", "description": "Name of the cell to stop"},
                    "title": {"type": "string", "description": "Title of the cell to stop"},
                    "action_summary": {"type": "string", "description": "Summary of the purpose of the stop."},
                },
                "required": ["cell_id", "title", "action_summary"],
            },
        })
        self.tool_map["stop_cell"] = stop_cell

        self.tools.append({
            "name": "delete_all_cells",
            "description": "Delete all cells in the notebook efficiently.",
            "input_schema": {
                "type": "object",
                "properties": {},
            },
        })
        self.tool_map["delete_all_cells"] = delete_all_cells

        self.tools.append({
            "name": "rename_notebook",
            "description": "Rename the current plot notebook.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "name": {"type": "string", "description": "New notebook name (<=5 words, Title Case)"},
                },
                "required": ["name"],
            },
        })
        self.tool_map["rename_notebook"] = rename_notebook

        self.tools.append({
            "name": "create_tab",
            "description": "Create a new tab at specified position to organize cells.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "position": {"type": "integer", "description": "Position to insert the tab marker"},
                    "display_name": {"type": "string", "description": "Name for the tab"},
                },
                "required": ["position", "display_name"],
            },
        })
        self.tool_map["create_tab"] = create_tab

        self.tools.append({
            "name": "rename_tab",
            "description": 'Rename a tab. Use tab_id="DEFAULT" to rename the default tab.',
            "input_schema": {
                "type": "object",
                "properties": {
                    "tab_id": {"type": "string", "description": 'ID of the tab to rename. Use "DEFAULT" for the default tab, or the TAB_ID from a Tab Marker.'},
                    "new_name": {"type": "string", "description": "New name for the tab"},
                },
                "required": ["tab_id", "new_name"],
            },
        })
        self.tool_map["rename_tab"] = rename_tab

        self.tools.append({
            "name": "submit_response",
            "description": "Submit the final response with plan, plan_diff, next_status, questions, and an optional summary. Call this at the end of every turn.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "plan": {"type": "array", "description": "List of plan items"},
                    "plan_diff": {"type": "array", "description": "List of plan diff items"},
                    "plan_update_overview": {"type": "string", "description": "Short title overview of what changed in the plan. Should follow the format like 'Added a new step.' or  `Completed step 2, step 3 now in progress.`"},
                    "summary": {"type": "string", "description": "Summary text to help the user. This can be a message to the user or a description of what was accomplished. Use markdown formatting with bullet points if needed. Omit if no summary needed."},
                    "questions": {"type": "string", "description": "Optional question text for the user. Omit if no questions needed."},
                    "next_status": {"type": "string", "description": "What the agent will do next", "enum": ["executing", "fixing", "thinking", "awaiting_user_response", "awaiting_cell_execution", "awaiting_user_widget_input", "done"]},
                    "expected_widgets": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Optional list of full widget keys (<tf_id>/<widget_id>) to await when next_status is 'awaiting_user_widget_input'"
                    },
                    "continue": {
                        "type": "boolean",
                        "description": "Set to true to immediately continue to the next step without waiting for user input. Set to false when waiting for user input or when all work is complete.",
                        "default": False
                    },
                },
                "required": ["plan", "plan_diff", "next_status", "plan_update_overview"],
            },
        })
        self.tool_map["submit_response"] = submit_response

        self.tools.append({
                    "name": "set_widget",
                    "description": "Set a single widget value by widget key.",
                    "input_schema": {
                        "type": "object",
                        "properties": {
                            "key": {
                                "type": "string",
                                "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>"
                            },
                            "value": {
                                "description": "JSON-serializable value"
                            },
                            "action_summary": {"type": "string", "description": "Summary of the purpose of the set_widget."},
                            "label": {"type": "string", "description": "Label of the widget to set"},
                        },
                        "required": ["key", "value", "action_summary", "label"],
                    },
                })
        self.tool_map["set_widget"] = set_widget

        self.tools.append({
            "name": "h5_filter_by",
            "description": "Set filters for an h5/AnnData widget. Pass the complete array of filters. Include existing filters from widget context to preserve them.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "widget_key": {
                        "type": "string",
                        "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>"
                    },
                    "label": {
                        "type": "string",
                        "description": "Label of the widget"
                    },
                    "filters": {
                        "type": "array",
                        "description": "Complete array of filters to apply",
                        "items": {
                            "oneOf": [
                                {
                                    "type": "object",
                                    "properties": {
                                        "type": {
                                            "type": "string",
                                            "enum": ["obs"],
                                            "description": "Filter by observation metadata"
                                        },
                                        "key": {
                                            "type": "string",
                                            "description": "The observation key to filter on"
                                        },
                                        "operation": {
                                            "oneOf": [
                                                {
                                                    "type": "object",
                                                    "properties": {
                                                        "type": {
                                                            "type": "string",
                                                            "enum": ["neq"],
                                                            "description": "Not equal operation"
                                                        },
                                                        "value": {
                                                            "type": ["string", "number", "null"],
                                                            "description": "Value to compare against"
                                                        }
                                                    },
                                                    "required": ["type", "value"]
                                                },
                                                {
                                                    "type": "object",
                                                    "properties": {
                                                        "type": {
                                                            "type": "string",
                                                            "enum": ["geq", "leq", "g", "l"],
                                                            "description": "Numeric comparison: geq (>=), leq (<=), g (>), l (<)"
                                                        },
                                                        "value": {
                                                            "type": "number",
                                                            "description": "Numeric value to compare against"
                                                        }
                                                    },
                                                    "required": ["type", "value"]
                                                }
                                            ],
                                            "description": "Filter operation to apply"
                                        }
                                    },
                                    "required": ["type", "key", "operation"]
                                },
                                {
                                    "type": "object",
                                    "properties": {
                                        "type": {
                                            "type": "string",
                                            "enum": ["var"],
                                            "description": "Filter by variable(s) / gene(s)"
                                        },
                                        "keys": {
                                            "type": "array",
                                            "items": {"type": "string"},
                                            "description": "Array of variable/gene names to filter on"
                                        },
                                        "operation": {
                                            "type": "object",
                                            "properties": {
                                                "type": {
                                                    "type": "string",
                                                    "enum": ["geq", "leq", "g", "l"],
                                                    "description": "Numeric comparison: geq (>=), leq (<=), g (>), l (<)"
                                                },
                                                "value": {
                                                    "type": "number",
                                                    "description": "Numeric value to compare against"
                                                }
                                            },
                                            "required": ["type", "value"]
                                        }
                                    },
                                    "required": ["type", "keys", "operation"]
                                }
                            ]
                        }
                    }
                },
                "required": ["widget_key", "filters"],
            },
        })
        self.tool_map["h5_filter_by"] = h5_filter_by

        self.tools.append({
            "name": "h5_color_by",
            "description": "Set an h5/AnnData widget to color by a specific observation or variable (can be multiple if for genes)",
            "input_schema": {
                "type": "object",
                "properties": {
                    "widget_key": {
                        "type": "string",
                        "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>"
                    },
                    "label": {
                        "type": "string",
                        "description": "Label of the widget"
                    },
                    "color_by": {
                        "oneOf": [
                            {
                                "type": "object",
                                "properties": {
                                    "type": {
                                        "type": "string",
                                        "enum": ["obs"],
                                        "description": "Color by observation metadata"
                                    },
                                    "key": {
                                        "type": "string",
                                        "description": "The observation key to color by"
                                    }
                                },
                                "required": ["type", "key"]
                            },
                            {
                                "type": "object",
                                "properties": {
                                    "type": {
                                        "type": "string",
                                        "enum": ["var"],
                                        "description": "Color by variable(s) / gene(s)"
                                    },
                                    "keys": {
                                        "type": "array",
                                        "items": {"type": "string"},
                                        "description": "Array of variable/gene names to color by"
                                    }
                                },
                                "required": ["type", "keys"]
                            },
                            {
                                "type": "null"
                            }
                        ],
                        "description": "Coloring configuration. Can be null to remove coloring, an obs object to color by observation, or a var object to color by variables like genes"
                    }
                },
                "required": ["widget_key", "color_by"],
            },
        })
        self.tool_map["h5_color_by"] = h5_color_by

        self.tools.append({
            "name": "h5_set_selected_obsm_key",
            "description": "Set the selected obsm key for an h5/AnnData widget to control which embedding is displayed.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "widget_key": {
                        "type": "string",
                        "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>"
                    },
                    "label": {
                        "type": "string",
                        "description": "Label of the widget"
                    },
                    "obsm_key": {
                        "type": "string",
                        "description": "The obsm key to use for embedding (e.g spatial, X_umap)"
                    },
                },
                "required": ["widget_key", "obsm_key"],
            },
        })
        self.tool_map["h5_set_selected_obsm_key"] = h5_set_selected_obsm_key

        self.tools.append({
            "name": "h5_set_background_image",
            "description": "Set a background image for an h5/AnnData widget from a user-attached file. The file must be an image type (jpg, jpeg, png, tiff).",
            "input_schema": {
                "type": "object",
                "properties": {
                    "label": {
                        "type": "string",
                        "description": "Label of the widget"
                    },
                    "widget_key": {
                        "type": "string",
                        "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>"
                    },
                    "node_id": {
                        "type": "string",
                        "description": "The LData node ID of the image file to use as background. This should come from files attached by the user in the chat."
                    },
                },
                "required": ["widget_key", "node_id"],
            },
        })
        self.tool_map["h5_set_background_image"] = h5_set_background_image

        self.tools.append({
            "name": "h5_open_image_aligner",
            "description": "Open the image alignment modal for a background image in an h5/AnnData widget.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "label": {
                        "type": "string",
                        "description": "Label of the widget"
                    },
                    "widget_key": {
                        "type": "string",
                        "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>"
                    },
                    "background_image_id": {
                        "type": "string",
                        "description": "The ID of the background image to align (typically the node_id)"
                    },
                },
                "required": ["widget_key", "background_image_id"],
            },
        })
        self.tool_map["h5_open_image_aligner"] = h5_open_image_aligner

        self.tools.append({
            "name": "h5_autoscale",
            "description": "Reset the plotted view in an h5/AnnData widget to Plotly's autoscaled data bounds.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "label": {
                        "type": "string",
                        "description": "Label of the widget"
                    },
                    "widget_key": {
                        "type": "string",
                        "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>"
                    },
                },
                "required": ["widget_key"],
            },
        })
        self.tool_map["h5_autoscale"] = h5_autoscale

        self.tools.append({
            "name": "h5_zoom",
            "description": "Zoom the Plotly view in an h5/AnnData widget in or out from the current camera center.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "widget_key": {
                        "type": "string",
                        "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>"
                    },
                    "label": {
                        "type": "string",
                        "description": "Label of the widget"
                    },
                    "direction": {
                        "type": "string",
                        "enum": ["in", "out"],
                        "description": 'Zoom direction; use "in" to zoom closer, "out" to zoom farther'
                    },
                    "percentage": {
                        "type": "number",
                        "minimum": 0,
                        "description": "Optional percentage change (e.g. 25 for ±25%); omitting uses the default Plotly zoom factor"
                    },
                },
                "required": ["widget_key", "direction"],
            },
        })
        self.tool_map["h5_zoom"] = h5_zoom

        self.tools.append({
            "name": "h5_set_background_image_visibility",
            "description": "Show or hide a specific background image in an h5/AnnData widget.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "label": {
                        "type": "string",
                        "description": "Label of the widget"
                    },
                    "widget_key": {
                        "type": "string",
                        "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>"
                    },
                    "background_image_id": {
                        "type": "string",
                        "description": "The ID of the background image to show or hide"
                    },
                    "hidden": {
                        "type": "boolean",
                        "description": "Whether to hide (true) or show (false) the background image"
                    },
                },
                "required": ["widget_key", "background_image_id", "hidden"],
            },
        })
        self.tool_map["h5_set_background_image_visibility"] = h5_set_background_image_visibility

        self.tools.append({
            "name": "h5_add_selected_cells_to_categorical_obs",
            "description": "Assign selected cells to a category in a categorical observation key",
            "input_schema": {
                "type": "object",
                "properties": {
                    "label": {
                        "type": "string",
                        "description": "Label of the widget"
                    },
                    "widget_key": {
                        "type": "string",
                        "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>"
                    },
                    "obs_key": {
                        "type": "string",
                        "description": "The existing categorical observation key to add selected cells to"
                    },
                    "category": {
                        "type": "string",
                        "description": "The category name to assign to selected cells. Will be created if it doesn't exist in this observation key."
                    },
                },
                "required": ["widget_key", "obs_key", "category"],
            },
        })
        self.tool_map["h5_add_selected_cells_to_categorical_obs"] = h5_add_selected_cells_to_categorical_obs

        self.tools.append({
            "name": "h5_set_marker_opacity",
            "description": "Set the marker opacity for all cell markers in an h5/AnnData widget.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "label": {
                        "type": "string",
                        "description": "Label of the widget"
                    },
                    "widget_key": {
                        "type": "string",
                        "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>"
                    },
                    "opacity": {
                        "type": "number",
                        "description": "Opacity value for cell markers, between 0.1 (transparent) and 0.9 (opaque)"
                    },
                },
                "required": ["widget_key", "opacity"],
            },
        })
        self.tool_map["h5_set_marker_opacity"] = h5_set_marker_opacity

        self.tools.append({
            "name": "h5_manage_obs",
            "description": "Create or delete an observation column in an h5/AnnData widget.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "widget_key": {
                        "type": "string",
                        "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>"
                    },
                    "label": {
                        "type": "string",
                        "description": "Label of the widget"
                    },
                    "obs_key": {
                        "type": "string",
                        "description": "The observation column name to create or delete"
                    },
                    "operation": {
                        "type": "string",
                        "enum": ["add", "remove"],
                        "description": "Whether to create a new observation column ('add') or delete an existing one ('remove')"
                    },
                    "obs_type": {
                        "type": "string",
                        "enum": ["category", "bool", "int64", "float64"],
                        "description": "Type of observation. Only for 'add' operation. Defaults to 'category'."
                    }
                },
                "required": ["widget_key", "obs_key", "operation"]
            }
        })
        self.tool_map["h5_manage_obs"] = h5_manage_obs

        self.tools.append({
            "name": "smart_ui_spotlight",
            "description": "Highlight a UI element to guide the user's attention.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "keyword": {
                        "type": "string",
                        "enum": ["lasso_select", "file_upload", "widget_input"],
                        "description": "The UI element to highlight"
                    },
                    "widget_key": {
                        "type": "string",
                        "description": "Optional full widget key including tf_id and widget_id in the format <tf_id>/<widget_id> for keywords related to a specific widget"
                    }
                },
                "required": ["keyword"]
            }
        })
        self.tool_map["smart_ui_spotlight"] = smart_ui_spotlight

        def glob_file_search(args: dict) -> dict:
            pattern = args.get("pattern", "")
            base_path = args.get("base_path", "agent_config/context")

            try:
                if not Path(base_path).is_absolute():
                    base_path = Path(__file__).parent / base_path
                else:
                    base_path = Path(base_path)

                if "/" in pattern:
                    flag = "-path"
                    search_pattern = f"*/{pattern}"
                else:
                    flag = "-name"
                    search_pattern = pattern

                result = subprocess.run(
                    ["/usr/bin/find", str(base_path), flag, search_pattern],
                    capture_output=True,
                    text=True,
                    timeout=5,
                    check=True
                )

                files = result.stdout.strip().split("\n") if result.stdout.strip() else []
                relative_files = [str(Path(f).relative_to(context_root)) for f in files if f]

                return {
                    "tool_name": "glob_file_search",
                    "success": True,
                    "summary": f"Found {len(relative_files)} files matching pattern '{pattern}'",
                    "files": relative_files,
                    "pattern": pattern,
                }
            except Exception as e:
                return {
                    "tool_name": "glob_file_search",
                    "success": False,
                    "summary": f"Error searching for files: {e}"
                }

        def grep(args: dict) -> dict:
            pattern = args.get("pattern")
            if pattern is None:
                return {
                    "tool_name": "grep",
                    "success": False,
                    "summary": "Pattern is required"
                }

            path = args.get("path", "agent_config/context")
            case_insensitive = args.get("case_insensitive", False)

            try:
                if not Path(path).is_absolute():
                    search_path = Path(__file__).parent / path
                else:
                    search_path = Path(path)

                cmd = ["/usr/bin/rg", "--line-number"]
                if case_insensitive:
                    cmd.append("--ignore-case")
                cmd.extend([pattern, str(search_path)])

                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=10,
                    check=False
                )

                display_path = str(search_path.relative_to(context_root))

                matches = result.stdout.strip().replace(str(context_root) + "/", "")

                if result.returncode == 0:
                    return {
                        "tool_name": "grep",
                        "success": True,
                        "summary": f"Found matches for pattern '{pattern}' in {display_path}",
                        "matches": matches,
                        "pattern": pattern,
                        "path": display_path
                    }

                if result.returncode == 1:
                    return {
                        "tool_name": "grep",
                        "success": True,
                        "summary": f"No matches found for pattern '{pattern}' in {display_path}",
                        "matches": "",
                        "pattern": pattern,
                        "path": display_path
                    }

                return {
                    "tool_name": "grep",
                    "success": False,
                    "summary": f"Ripgrep error: {result.stderr}"
                }
            except Exception as e:
                return {
                    "tool_name": "grep",
                    "success": False,
                    "summary": f"Error searching: {e}"
                }

        def read_file(args: dict) -> dict:
            path = args.get("path")
            if path is None:
                return {
                    "tool_name": "read_file",
                    "success": False,
                    "summary": "Path is required"
                }

            offset = args.get("offset", 0)
            limit = args.get("limit")

            try:
                if not Path(path).is_absolute():
                    file_path = Path(__file__).parent / path
                else:
                    file_path = Path(path)

                if not file_path.exists():
                    return {
                        "tool_name": "read_file",
                        "success": False,
                        "summary": f"File not found: {path}"
                    }

                with Path(file_path).open(encoding="utf-8") as f:
                    lines = f.readlines()

                total_lines = len(lines)

                if offset >= total_lines:
                    return {
                        "tool_name": "read_file",
                        "success": False,
                        "summary": f"Offset {offset} exceeds file length {total_lines}"
                    }

                if limit is not None:
                    selected_lines = lines[offset:offset + limit]
                else:
                    selected_lines = lines[offset:]

                numbered_lines = []
                for i, line in enumerate(selected_lines, start=offset + 1):
                    numbered_lines.append(f"{i:6}|{line.rstrip()}")

                content = "\n".join(numbered_lines)

                display_path = str(file_path.relative_to(context_root))

                return {
                    "tool_name": "read_file",
                    "success": True,
                    "summary": f"Read {len(selected_lines)} lines from {display_path} (total: {total_lines} lines)",
                    "content": content,
                    "path": display_path,
                    "offset": offset,
                    "lines_read": len(selected_lines),
                    "total_lines": total_lines
                }
            except Exception as e:
                return {
                    "tool_name": "read_file",
                    "success": False,
                    "summary": f"Error reading file: {e}"
                }

        def search_replace(args: dict) -> dict:
            path = args.get("path")
            old_string = args.get("old_string")
            new_string = args.get("new_string", "")

            if path is None or old_string is None:
                return {
                    "tool_name": "search_replace",
                    "success": False,
                    "summary": "Path and old_string are required"
                }

            try:
                if not Path(path).is_absolute():
                    file_path = Path(__file__).parent / path
                else:
                    file_path = Path(path)

                if not file_path.exists():
                    return {
                        "tool_name": "search_replace",
                        "success": False,
                        "summary": f"File not found: {path}"
                    }

                tmp_file = Path(str(file_path) + ".tmp")

                result = subprocess.run(  # noqa: S602
                    f"rg --passthru --fixed-strings '{old_string}' --replace '{new_string}' {file_path} > {tmp_file} && mv {tmp_file} {file_path}",
                    shell=True,
                    capture_output=True,
                    text=True,
                    check=False
                )

                display_path = str(file_path.relative_to(context_root))

                if result.returncode == 0:
                    return {
                        "tool_name": "search_replace",
                        "success": True,
                        "summary": f"Replaced text in {display_path}",
                        "path": display_path
                    }

                if tmp_file.exists():
                    tmp_file.unlink()

                return {
                    "tool_name": "search_replace",
                    "success": False,
                    "summary": f"Ripgrep replace failed: {result.stderr}"
                }
            except Exception as e:
                return {
                    "tool_name": "search_replace",
                    "success": False,
                    "summary": f"Error replacing text: {e}"
                }

        def bash(args: dict) -> dict:
            command = args.get("command", "")

            try:
                result = subprocess.run(  # noqa: S602
                    command,
                    check=False,
                    shell=True,
                    capture_output=True,
                    text=True,
                    timeout=30,
                    cwd=str(Path(__file__).parent / "agent_config/context")
                )

                output = result.stdout + result.stderr

                return {
                    "tool_name": "bash",
                    "success": result.returncode == 0,
                    "summary": f"Executed command: {command}",
                    "output": output.strip(),
                    "command": command,
                    "return_code": result.returncode
                }
            except subprocess.TimeoutExpired:
                return {
                    "tool_name": "bash",
                    "success": False,
                    "summary": "Command timed out after 30 seconds"
                }
            except Exception as e:
                return {
                    "tool_name": "bash",
                    "success": False,
                    "summary": f"Error executing command: {e}"
                }

        self.tools.append({
            "name": "glob_file_search",
            "description": "Search for files matching a pattern. Defaults to searching the context directory. For nested paths use patterns like 'technology_docs/*.md'.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "pattern": {
                        "type": "string",
                        "description": "Glob pattern to match files. Use '*.md' for all markdown files, 'technology_docs/*.md' for files in a subdirectory."
                    },
                    "base_path": {
                        "type": "string",
                        "description": "Base path to search in. Defaults to 'agent_config/context'. Use relative paths from agent directory."
                    }
                },
                "required": ["pattern"]
            }
        })
        self.tool_map["glob_file_search"] = glob_file_search

        self.tools.append({
            "name": "grep",
            "description": "Search for text patterns in files using ripgrep (rg). Returns matches with line numbers. Fast and supports regex.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "pattern": {
                        "type": "string",
                        "description": "Text pattern to search for (supports regex)"
                    },
                    "path": {
                        "type": "string",
                        "description": "File or directory path to search in. Defaults to 'agent_config/context'."
                    },
                    "case_insensitive": {
                        "type": "boolean",
                        "description": "Whether to perform case-insensitive search. Defaults to false."
                    }
                },
                "required": ["pattern"]
            }
        })
        self.tool_map["grep"] = grep

        self.tools.append({
            "name": "read_file",
            "description": "Read contents of a file with optional offset and limit for large files. Returns numbered lines.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Path to the file to read (relative to agent directory or absolute)"
                    },
                    "offset": {
                        "type": "integer",
                        "description": "Line number to start reading from (0-indexed). Defaults to 0."
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of lines to read. If not provided, reads to end of file."
                    }
                },
                "required": ["path"]
            }
        })
        self.tool_map["read_file"] = read_file

        self.tools.append({
            "name": "search_replace",
            "description": "Replace the first occurrence of a string in a file with another string. Useful for editing files or maintaining state.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Path to the file to edit"
                    },
                    "old_string": {
                        "type": "string",
                        "description": "Exact string to find and replace"
                    },
                    "new_string": {
                        "type": "string",
                        "description": "String to replace with"
                    }
                },
                "required": ["path", "old_string", "new_string"]
            }
        })
        self.tool_map["search_replace"] = search_replace

        self.tools.append({
            "name": "bash",
            "description": "Execute a bash command. The working directory is already set to the context directory, so use relative paths (e.g., 'ls technology_docs/' not 'ls agent_config/context/technology_docs/').",
            "input_schema": {
                "type": "object",
                "properties": {
                    "command": {
                        "type": "string",
                        "description": "Bash command to execute. Working directory is the context directory, use relative paths."
                    }
                },
                "required": ["command"]
            }
        })
        self.tool_map["bash"] = bash

        async def execute_code(args: dict) -> dict:
            code = args.get("code")
            if code is None:
                return {
                    "tool_name": "execute_code",
                    "success": False,
                    "summary": "No code provided"
                }

            print(f"[tool] execute_code: {code[:50]}...")

            result = await self.atomic_operation("execute_code", {"code": code})

            return {
                "tool_name": "execute_code",
                "success": True,
                "summary": "Code executed",
                "code": code,
                "stdout": result.get("stdout"),
                "stderr": result.get("stderr"),
                "exception": result.get("exception"),
            }

        self.tools.append({
            "name": "execute_code",
            "description": "Execute arbitrary Python code in the notebook kernel and return the result, stdout, stderr, and any exceptions. Use this to test imports, print values, or run simple inspection code.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "code": {"type": "string", "description": "Python code to execute"}
                },
                "required": ["code"]
            }
        })
        self.tool_map["execute_code"] = execute_code

        async def get_global_info(args: dict) -> dict:
            key = args.get("key")
            if key is None:
                return {
                    "tool_name": "get_global_info",
                    "success": False,
                    "summary": "No key provided"
                }

            print(f"[tool] get_global_info: {key}")

            result = await self.atomic_operation("get_global_info", {"key": key})

            if result.get("status") == "success":
                info = result.get("info", {})
                return {
                    "tool_name": "get_global_info",
                    "success": True,
                    "summary": f"Retrieved info for global '{key}'",
                    "key": key,
                    "info": info
                }

            return {
                "tool_name": "get_global_info",
                "success": False,
                "summary": f"Failed to get global info: {result.get('error', 'Unknown error')}"
            }

        self.tools.append({
            "name": "get_global_info",
            "description": "Get rich information about a specific global variable including its type, shape, columns, dtypes, etc. Especially useful for DataFrames and AnnData objects.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "key": {"type": "string", "description": "Name of the global variable to inspect"}
                },
                "required": ["key"]
            }
        })
        self.tool_map["get_global_info"] = get_global_info

        async def refresh_cells_context(args: dict) -> dict:
            context_result = await self.atomic_operation("get_context")

            if context_result.get("status") != "success":
                return {
                    "tool_name": "refresh_cells_context",
                    "success": False,
                    "summary": f"Failed to refresh cells: {context_result.get('error', 'Unknown error')}"
                }

            context = context_result.get("context", {})
            self.latest_notebook_context = context

            notebook_name = context.get("notebook_name")
            cell_count = context.get("cell_count", 0)
            cells = context.get("cells", [])

            cell_lines = [f"# Notebook Cells for {notebook_name}, Total cells: {cell_count}\n"]

            default_tab_name = context.get("default_tab_name", "Tab 1")

            cell_lines.append("\n## Tab Marker [DEFAULT]")  # noqa: FURB113
            cell_lines.append(f"TAB_NAME: {default_tab_name}")
            cell_lines.append("TAB_ID: DEFAULT")
            cell_lines.append("TYPE: Default Tab Marker")
            cell_lines.append("---")

            current_tab_name = default_tab_name

            for cell in cells:
                index = cell.get("index", "?")
                cell_id = cell.get("cell_id", "?")
                cell_type = cell.get("cell_type", "unknown")
                source = cell.get("source")
                status = cell.get("status", "idle")
                tf_id = cell.get("tf_id", None)

                if cell_type == "tabMarker":
                    if source:
                        current_tab_name = source.strip() or f"Tab {index}"
                    else:
                        current_tab_name = f"Tab {index}"

                    cell_lines.append(f"\n## Tab Marker [{index}]")  # noqa: FURB113
                    cell_lines.append(f"TAB_NAME: {current_tab_name}")
                    cell_lines.append(f"TAB_ID: {cell_id}")
                    cell_lines.append(f"CELL_INDEX: {index}")
                    cell_lines.append("TYPE: Tab Marker")
                    cell_lines.append("---")
                    continue

                cell_lines.append(f"\n## Cell [{index}] (in {current_tab_name})")  # noqa: FURB113
                cell_lines.append(f"BELONGS_TO_TAB: {current_tab_name}")
                cell_lines.append(f"CELL_ID: {cell_id}")
                cell_lines.append(f"CELL_INDEX: {index}")
                cell_lines.append(f"TYPE: {cell_type}")
                cell_lines.append(f"STATUS: {status}")
                if tf_id is not None:
                    cell_lines.append(f"CODE_CELL_ID: {tf_id}")

                if source is not None:
                    cell_lines.append("CODE_START")  # noqa: FURB113
                    cell_lines.append(source)
                    cell_lines.append("CODE_END")

                widgets = cell.get("widgets")
                if widgets is not None:
                    cell_lines.append("\nWIDGETS:")
                    for w in widgets:
                        w_type = w.get("type", "unknown")
                        w_key = w.get("key", "")
                        w_label = w.get("label", "")
                        cell_lines.append(f"- WIDGET: {w_type} | {w_label} | {w_key}")

            context_dir = context_root / "notebook_context"
            context_dir.mkdir(parents=True, exist_ok=True)
            (context_dir / "cells.md").write_text("\n".join(cell_lines))

            return {
                "tool_name": "refresh_cells_context",
                "success": True,
                "summary": f"Refreshed cells context for {cell_count} cells and stored result in {context_dir / 'cells.md'}",
                "cell_count": cell_count,
                "context_path": str(context_dir / "cells.md")
            }

        async def refresh_reactivity_context(args: dict) -> dict:
            reactivity_result = await self.atomic_operation("request_reactivity_summary")

            if reactivity_result.get("status") != "success":
                return {
                    "tool_name": "refresh_reactivity_context",
                    "success": False,
                    "summary": f"Failed to refresh reactivity: {reactivity_result.get('error', 'Unknown error')}"
                }

            reactivity_summary = reactivity_result.get("summary")

            context_dir = context_root / "notebook_context"
            context_dir.mkdir(parents=True, exist_ok=True)

            if reactivity_summary is not None:
                (context_dir / "signals.md").write_text(reactivity_summary)
            else:
                (context_dir / "signals.md").write_text("# Reactivity Summary\n\nNo reactive dependencies in this notebook.\n")

            return {
                "tool_name": "refresh_reactivity_context",
                "success": True,
                "summary": f"Refreshed reactivity context and stored result in {context_dir / 'signals.md'}",
                "context_path": str(context_dir / "signals.md"),
            }

        self.tools.append({
            "name": "refresh_cells_context",
            "description": "Refresh the cells.md context file with current notebook cell structure and contents.",
            "input_schema": {
                "type": "object",
                "properties": {},
            },
        })
        self.tool_map["refresh_cells_context"] = refresh_cells_context

        self.tools.append({
            "name": "refresh_reactivity_context",
            "description": "Refresh the signals.md context file with current reactive signal dependencies.",
            "input_schema": {
                "type": "object",
                "properties": {},
            },
        })
        self.tool_map["refresh_reactivity_context"] = refresh_reactivity_context

        if len(self.tools) > 0:
            self.tools[-1]["cache_control"] = {"type": "ephemeral"}

    async def handle_stream(self, *, use_beta_api: bool, **kwargs) -> tuple[Message | BetaMessage, float]:  # noqa: ANN003
        start_time = time.time()

        await self.send({
            "type": "agent_stream_start",
            "timestamp": int(time.time() * 1000),
        })

        content_blocks: list[dict] = []
        current_block_index = -1
        usage_data = None

        try:
            stream_ctx = self.client.beta.messages.stream(**kwargs) if use_beta_api else self.client.messages.stream(**kwargs)

            async with stream_ctx as stream:
                async for event in stream:
                    event_type = event.type

                    if event_type == "message_start":
                        message_data = event.message
                        if hasattr(message_data, "usage") and message_data.usage:
                            usage_data = message_data.usage

                    elif event_type == "content_block_start":
                        current_block_index = event.index
                        block = event.content_block

                        block_dict = {"type": block.type}
                        if block.type == "text":
                            block_dict["text"] = ""
                        elif block.type == "thinking":
                            block_dict["thinking"] = ""
                        elif block.type == "tool_use":
                            block_dict["id"] = block.id
                            block_dict["name"] = block.name
                            block_dict["input"] = {}

                        content_blocks.append(block_dict)

                        block_start_msg = {
                            "type": "agent_stream_block_start",
                            "block_index": current_block_index,
                            "block_type": block.type,
                        }
                        if block.type == "tool_use":
                            block_start_msg["block_id"] = block.id
                            block_start_msg["block_name"] = block.name

                        await self.send(block_start_msg)

                    elif event_type == "content_block_delta":
                        delta = event.delta
                        block_index = event.index

                        if delta.type == "text_delta":
                            content_blocks[block_index]["text"] += delta.text
                            await self.send({
                                "type": "agent_stream_delta",
                                "block_index": block_index,
                                "block_type": "text",
                                "delta": delta.text,
                            })

                        elif delta.type == "thinking_delta":
                            content_blocks[block_index]["thinking"] += delta.thinking
                            await self.send({
                                "type": "agent_stream_delta",
                                "block_index": block_index,
                                "block_type": "thinking",
                                "delta": delta.thinking,
                            })

                        elif delta.type == "input_json_delta":
                            partial_json = delta.partial_json
                            await self.send({
                                "type": "agent_stream_delta",
                                "block_index": block_index,
                                "block_type": "tool_use",
                                "delta": partial_json,
                            })

                    elif event_type == "content_block_stop":
                        block_index = event.index
                        await self.send({
                            "type": "agent_stream_block_stop",
                            "block_index": block_index,
                        })

                    elif event_type == "message_delta":
                        if hasattr(event, "usage") and event.usage:
                            usage_data = event.usage

                    elif event_type == "message_stop":
                        pass

                final_message = await stream.get_final_message()

            duration_seconds = time.time() - start_time

            if usage_data is not None:
                cache_read_input_tokens = getattr(usage_data, "cache_read_input_tokens", None) or 0
                cache_creation_input_tokens = getattr(usage_data, "cache_creation_input_tokens", None) or 0

                await self.send({
                    "type": "agent_usage_update",
                    "input_tokens": usage_data.input_tokens,
                    "cache_read_input_tokens": cache_read_input_tokens,
                    "cache_creation_input_tokens": cache_creation_input_tokens,
                    "context_limit": 200_000,  # todo(aidan): store this info in db per model config
                })

            await self.send({
                "type": "agent_stream_complete",
            })

            print(f"[agent] Stream completed in {duration_seconds:.3f}s, {len(content_blocks)} blocks")

            return final_message, duration_seconds

        except Exception as e:
            print(f"[agent] Stream error: {e}")
            traceback.print_exc()

            await self.send({
                "type": "agent_stream_complete",
                "error": str(e),
            })

            raise

    async def run_agent_loop(self) -> None:
        assert self.client is not None, "Client not initialized"

        self.conversation_running = True
        turn = 0

        print("[agent] run_agent_loop: started")

        while self.conversation_running:
            await self._wait_for_message()
            if not self.conversation_running:
                break

            print("[agent] run_agent_loop: building messages from DB...")
            build_start = time.time()
            api_messages = self._truncate_old_messages(await self._build_messages_from_db())
            build_elapsed = time.time() - build_start
            print(f"[agent] run_agent_loop: built {len(api_messages) if api_messages else 0} messages in {build_elapsed:.3f}s")

            if not api_messages or api_messages[-1].get("role") != "user":
                print("[agent] run_agent_loop: skipping (no messages or last message not user)")
                continue

            turn += 1
            print(f"[agent] run_agent_loop: starting turn {turn}")

            model, thinking_budget = self.mode_config.get(self.mode, ("claude-sonnet-4-5-20250929", 1024))

            print(f"[agent] Turn {turn}, mode={self.mode}, thinking_budget={thinking_budget}")

            if thinking_budget is not None:
                max_tokens = thinking_budget + 4096
            else:
                max_tokens = 4096

            can_use_thinking = True
            if thinking_budget is not None and len(api_messages) > 0:
                has_any_thinking = False
                for msg in api_messages:
                    if msg.get("role") == "assistant":
                        content = msg.get("content", [])
                        if isinstance(content, list):
                            for block in content:
                                block_type = block.get("type") if isinstance(block, dict) else None
                                if block_type in {"thinking", "redacted_thinking"}:
                                    has_any_thinking = True
                                    break
                        if has_any_thinking:
                            break

                if not has_any_thinking:
                    last_assistant_msg = None
                    for msg in reversed(api_messages):
                        if msg.get("role") == "assistant":
                            last_assistant_msg = msg
                            break

                    if last_assistant_msg:
                        content = last_assistant_msg.get("content", [])
                        if isinstance(content, list) and len(content) > 0:
                            first_block_type = content[0].get("type") if isinstance(content[0], dict) else None
                            if first_block_type not in {"thinking", "redacted_thinking"}:
                                can_use_thinking = False
                                print(f"[agent] Cannot use thinking API: last assistant message starts with {first_block_type}, not thinking")

            assert self.system_prompt is not None

            system_blocks = [
                {
                    "type": "text",
                    "text": self.system_prompt,
                    "cache_control": {"type": "ephemeral"}
                }
            ]

            kwargs = {
                "model": model,
                "max_tokens": max_tokens,
                "system": system_blocks,
                "messages": api_messages,
                "tools": self.tools,
            }
            use_beta_api = False
            if thinking_budget is not None and can_use_thinking:
                kwargs["thinking"] = {
                    "type": "enabled",
                    "budget_tokens": thinking_budget,
                }
                kwargs["betas"] = ["interleaved-thinking-2025-05-14"]
                use_beta_api = True

            try:
                response, duration_seconds = await self.handle_stream(
                    use_beta_api=use_beta_api,
                    **kwargs
                )

            except Exception as e:
                print(f"[agent] API error: {e}")

                print(f"[agent] API call failed with {len(api_messages)} messages")
                for i, msg in enumerate(api_messages):
                    role = msg.get("role", "?")
                    content = msg.get("content", [])
                    if isinstance(content, str):
                        print(f"  Message {i} ({role}): string content, length={len(content)}")
                    elif isinstance(content, list):
                        print(f"  Message {i} ({role}): {len(content)} blocks")
                        for j, block in enumerate(content):
                            block_type = block.get("type") if isinstance(block, dict) else getattr(block, "type", "?")
                            print(f"    Block {j}: {block_type}")
                    else:
                        print(f"  Message {i} ({role}): unknown content type={type(content)}")

                await self.send({
                    "type": "agent_error",
                    "error": f"API error: {e!s}",
                    "fatal": False
                })
                continue

            response_content = response.model_dump()["content"]
            if response_content is not None and (not isinstance(response_content, list) or len(response_content) > 0):
                await self._insert_history(
                    role="assistant",
                    payload={
                        "content": response_content,
                        "duration": duration_seconds,
                    },
                )
            else:
                print(f"[agent] Skipping empty assistant message (stop_reason={response.stop_reason})")

            if response.stop_reason == "end_turn":
                print("[agent] Turn ended without submit_response; completing turn")
                await self._complete_turn()
            elif response.stop_reason == "tool_use":
                tool_results = []
                called_submit_response = False

                for block in response.content:
                    if isinstance(block, dict):
                        block_type = block.get("type")
                    else:
                        block_type = getattr(block, "type", None)

                    if block_type == "tool_use":
                        tool_id = block.get("id") if isinstance(block, dict) else block.id
                        tool_name = block.get("name") if isinstance(block, dict) else block.name
                        tool_input = block.get("input") if isinstance(block, dict) else block.input

                        if tool_name == "submit_response":
                            called_submit_response = True

                        print(f"[agent] Executing tool: {tool_name} (id={tool_id})")

                        handler = self.tool_map.get(tool_name)

                        self.pending_tool_calls.add(tool_id)
                        if handler is not None:
                            try:
                                result = handler(tool_input)
                                if asyncio.iscoroutine(result):
                                    result = await result

                                tool_results.append({
                                    "type": "tool_result",
                                    "tool_use_id": tool_id,
                                    "content": json.dumps(result),
                                })
                            except Exception as e:
                                print(f"[agent] Tool error: {tool_name}: {e}")
                                traceback.print_exc()
                                tool_results.append({
                                    "type": "tool_result",
                                    "tool_use_id": tool_id,
                                    "content": json.dumps({
                                        "tool_name": tool_name,
                                        "summary": f"Error executing tool: {e!s}",
                                        "error": str(e),
                                        "success": False,
                                    }),
                                })
                        else:
                            tool_results.append({
                                "type": "tool_result",
                                "tool_use_id": tool_id,
                                "content": json.dumps({
                                    "summary": None,
                                    "error": f"Unknown tool: {tool_name}",
                                }),
                            })

                if len(tool_results) > 0:
                    uncancelled_tool_results = [result for result in tool_results if result.get("tool_use_id") in self.pending_tool_calls]
                    self.pending_tool_calls.difference_update({result.get("tool_use_id") for result in tool_results})
                    await self._insert_history(
                        payload={
                            "content": uncancelled_tool_results,
                        },
                    )
                else:
                    print("[agent] No tool results")

                if called_submit_response:
                    print("[agent] submit_response called, completing turn")
                    await self._complete_turn()
                else:
                    await self.pending_messages.put({"type": "resume"})
            elif response.stop_reason == "max_tokens":
                print("[agent] Hit max tokens")
                await self._complete_turn()
            else:
                print(f"[agent] Unknown stop reason: {response.stop_reason}")
                await self._complete_turn()

    async def _close_pending_tool_calls(self, *, error_message: str, messages: list[MessageParam]) -> None:
        tool_use_ids = set()
        tool_result_ids = set()

        for history_msg in messages:
            content = history_msg.get("content")
            if not isinstance(content, list):
                continue

            for block in content:
                if not isinstance(block, dict):
                    continue

                if block.get("type") == "tool_use":
                    tool_id = block.get("id")
                    if tool_id:
                        tool_use_ids.add(tool_id)
                elif block.get("type") == "tool_result":
                    tool_use_id = block.get("tool_use_id")
                    if tool_use_id:
                        tool_result_ids.add(tool_use_id)

        pending_tool_ids = tool_use_ids - tool_result_ids

        if len(pending_tool_ids) > 0:
            print(f"[agent] Closing {len(pending_tool_ids)} pending tool calls")
            await self._insert_history(
                payload={
                    "content": [
                        {
                            "type": "tool_result",
                            "tool_use_id": tool_id,
                            "content": json.dumps({
                                "summary": None,
                                "error": error_message,
                            }),
                        }
                        for tool_id in pending_tool_ids
                    ],
                },
            )
            await self._notify_history_updated()

    _context_init_task: asyncio.Task | None = None

    async def handle_init(self, msg: dict[str, object]) -> None:
        print("[agent] Initializing")

        try:
            await self._clear_running_state()

            self.should_auto_continue = False
            self.pending_auto_continue = False

            session_id = msg.get("session_id")
            if session_id is None:
                raise RuntimeError(f"[handle init] Session ID is not set. Message: {msg}")

            self.agent_session_id = int(session_id)

            self.init_tools()

            self.client = anthropic.AsyncAnthropic(
                api_key="dummy",
                base_url=f"{nucleus_url}/infer/plots-agent/anthropic",
                default_headers={"Authorization": auth_token_sdk, "Pod-Id": str(pod_id)}
            )

            system_prompt_path = context_root.parent / "system_prompt.md"
            self.system_prompt = system_prompt_path.read_text()

            self.pending_tool_calls.clear()
            messages = await self._build_messages_from_db()
            await self._close_pending_tool_calls(
                error_message="Tool call cancelled - session was ended before completion",
                messages=messages,
            )

            self.initialized = True
            await self.send({
                "type": "agent_status",
                "status": "ready"
            })
            print("[agent] Initialization complete")

            print("[agent] Starting conversation loop")
            self._start_conversation_loop()
            print("[agent] Conversation loop started")

            most_recent_submit_response = None
            for history_msg in reversed(messages):
                if history_msg.get("role") == "assistant":
                    content = history_msg.get("content")
                    if isinstance(content, list):
                        for block in content:
                            if isinstance(block, dict) and block.get("type") == "tool_use" and block.get("name") == "submit_response":
                                most_recent_submit_response = block.get("input", {})
                                break
                    if most_recent_submit_response is not None:
                        break

            if most_recent_submit_response is not None:
                next_status = most_recent_submit_response.get("next_status")
                waiting_states = {"awaiting_cell_execution", "awaiting_user_widget_input"}

                if next_status in waiting_states:
                    previous_request_id = await self._extract_last_request_id()
                    print(f"[agent] Reconnected while {next_status}, prompting LLM to check state and retry (request_id={previous_request_id})")
                    await self.pending_messages.put({
                        "type": "user_query",
                        "content": "The session was reconnected. You were waiting for an action to complete, but it may have finished or failed while offline. Please retry any incomplete actions or continue with your plan.",
                        "request_id": previous_request_id,
                        "hidden": True,
                    })
                else:
                    print(f"[agent] Reconnected with status '{next_status}', staying idle")

            if len(messages) > 0 and messages[-1].get("role") == "user":
                print("[agent] Incomplete turn detected, auto-resuming")
                await self.pending_messages.put({"type": "resume"})
        except Exception as e:
            await self.send({
                "type": "agent_error",
                "error": f"Failed to initialize: {e!s}",
                "fatal": True
            })

    async def handle_query(self, msg: dict[str, object]) -> None:
        query = msg.get("query", "")
        request_id = msg.get("request_id")
        contextual_node_data = msg.get("contextual_node_data")

        full_query = query
        if contextual_node_data:
            full_query = f"{query} \n\nHere is the context of the selected nodes the user would like to use: <ContextualNodeData>{json.dumps(contextual_node_data)}</ContextualNodeData>"

        await self.pending_messages.put({
            "type": "user_query",
            "content": full_query,
            "request_id": request_id,
            "display_query": query,
            "display_nodes": contextual_node_data,
        })

        print(f"[agent] Message queued successfully (request_id={request_id})")

    async def handle_cancel(self, msg: dict[str, object]) -> None:
        request_id = msg.get("request_id", "unknown")
        print(f"[agent] Cancelling request {request_id}")

        await self._clear_running_state()

        self.current_request_id = None
        self.should_auto_continue = False
        self.pending_auto_continue = False
        self.executing_cells.clear()

        self.pending_tool_calls.clear()
        await self._close_pending_tool_calls(
            error_message="Request cancelled by user",
            messages=await self._build_messages_from_db(),
        )

        for file in (context_root / "agent_scratch").rglob("*"):
            if file.name == ".gitkeep":
                continue

            file.unlink()

        self._start_conversation_loop()

    async def handle_clear_history(self) -> None:
        await self._clear_running_state()
        await self._mark_all_history_removed()
        self._start_conversation_loop()

    async def get_full_prompt(self) -> dict:
        messages = await self._build_messages_from_db()

        context_dir = context_root / "notebook_context"

        def read_context_file(filename: str) -> str:
            file_path = context_dir / filename
            if file_path.exists():
                return file_path.read_text()
            return f"# {filename}\n\nFile not yet generated."

        cells_content = read_context_file("cells.md")
        signals_content = read_context_file("signals.md")

        def build_tree(path: Path, prefix: str = "") -> list[str]:
            lines = []
            entries = sorted(path.iterdir(), key=lambda p: (not p.is_dir(), p.name))
            for entry in entries:
                if entry.is_dir():
                    lines.append(f"{prefix}{entry.name}/")
                    lines.extend(build_tree(entry, prefix + "  "))
                else:
                    lines.append(f"{prefix}{entry.name}")

            return lines

        tree_lines = ["# Context Directory Structure", "", "context/"]
        tree_lines.extend(build_tree(context_root, "  "))
        tree_content = "\n".join(tree_lines)

        truncated_messages = self._truncate_old_messages(messages)

        return {
            "system_prompt": self.system_prompt,
            "messages": messages,
            "truncated_messages": truncated_messages,
            "model": self.mode_config.get(self.mode, ("claude-sonnet-4-5-20250929", 1024))[0],
            "cells": cells_content,
            "signals": signals_content,
            "tree": tree_content,
        }

    async def update_system_prompt(self, msg: dict[str, object]) -> dict:
        new_content = msg.get("content")
        if not isinstance(new_content, str):
            return {"status": "error", "error": "Invalid content"}

        system_prompt_path = context_root.parent / "system_prompt.md"
        system_prompt_path.write_text(new_content)
        self.system_prompt = new_content

        full_prompt = await self.get_full_prompt()
        return {
            "status": "success",
            **full_prompt,
        }

    async def accept(self) -> None:
        msg = await self.conn.recv()
        msg_type = msg.get("type")
        msg_request_id = msg.get("request_id")

        print(f"[agent] accept: received message type={msg_type} (request_id={msg_request_id})")

        if msg_type == "init":
            print(f"[agent] Message: {msg_type}")
            await self.handle_init(msg)
        elif msg_type == "agent_query":
            query = msg.get("query", "")
            query_preview = query[:60] + "..." if len(query) > 60 else query
            request_id = msg.get("request_id", "unknown")
            print(f"[agent] accept: dispatching to handle_query (query={query_preview}, request_id={request_id})")
            handle_start = time.time()
            await self.handle_query(msg)
            handle_elapsed = time.time() - handle_start
            print(f"[agent] accept: handle_query completed in {handle_elapsed:.3f}s (request_id={request_id})")
        elif msg_type == "agent_cancel":
            request_id = msg.get("request_id", "unknown")
            print(f"[agent] Cancel: {request_id}")
            await self.handle_cancel(msg)
        elif msg_type == "agent_clear_history":
            print("[agent] Clear history request")
            await self.handle_clear_history()
        elif msg_type == "agent_action_response":
            print(f"[agent] {msg.get('action', 'unknown')} -> {msg.get('status', 'unknown')}")
            await self.handle_action_response(msg)
        elif msg_type == "kernel_message":
            print(f"[agent] Kernel message: {msg}")

            nested_msg = msg.get("message", {})
            nested_type = nested_msg.get("type")

            if nested_type == "cell_result":
                cell_id = nested_msg.get("cell_id")
                has_exception = nested_msg.get("has_exception", False)
                exception = nested_msg.get("exception")
                display_name = nested_msg.get("display_name")

                logs = nested_msg.get("logs", None)
                if logs is not None and len(logs) > 4096:
                    logs = logs[-4096:]

                if cell_id is not None:
                    self.executing_cells.discard(str(cell_id))
                if self.current_status == "awaiting_user_widget_input" and not all(v is None for v in self.expected_widgets.values()):
                    print(f"[agent] Not adding cell {cell_id} result because awaiting_user_widget_input and still expecting widget updates")
                    return

                if self.current_request_id is not None:
                    await self.pending_messages.put({
                        "type": "cell_result",
                        "cell_id": cell_id,
                        "success": not has_exception,
                        "exception": exception,
                        "logs": logs,
                        "display_name": display_name,
                    })
                else:
                    print(f"[agent] Cell {cell_id} completed but no active request - updating executing_cells only")

            elif nested_type == "start_cell":
                cell_id = nested_msg.get("cell_id")
                if cell_id is not None and self.current_request_id is not None:
                    self.executing_cells.add(str(cell_id))
            elif nested_type == "set_widget_value":
                if self.current_status == "awaiting_user_widget_input":
                    data = nested_msg.get("data", {})
                    for key, value in data.items():
                        if key in self.expected_widgets:
                            self.expected_widgets[key] = value

                    if all(v is not None for v in self.expected_widgets.values()):
                        await self.pending_messages.put({
                            "type": "set_widget_value",
                            "data": self.expected_widgets
                        })
        elif msg_type == "get_full_prompt":
            tx_id = msg.get("tx_id")
            print(f"[agent] Get full prompt request (tx_id={tx_id})")

            result = await self.get_full_prompt()
            await self.send({
                "type": "agent_action_response",
                "tx_id": tx_id,
                "status": "success",
                **result
            })
        elif msg_type == "update_system_prompt":
            tx_id = msg.get("tx_id")
            print(f"[agent] Update system prompt request (tx_id={tx_id})")
            result = await self.update_system_prompt(msg)

            await self.send({
                "type": "agent_action_response",
                "tx_id": tx_id,
                **result
            })
        else:
            print(f"[agent] Unknown message type: {msg_type}")


async def main() -> None:
    global loop
    loop = asyncio.get_running_loop()

    from datetime import datetime
    print(f"{datetime.now().isoformat()} [agent] Starting")

    sock = socket.socket(family=socket.AF_UNIX, fileno=int(sys.argv[-1]))
    sock.setblocking(False)

    socket_io_thread = SocketIoThread(socket=sock)
    socket_io_thread.start()
    try:
        socket_io_thread.initialized.wait()

        harness = AgentHarness(conn=socket_io_thread)
        _inject.agent = harness

        await harness.send({"type": "ready"})

        while True:
            try:
                await harness.accept()
            except Exception:
                traceback.print_exc()

        print("Agent shutting down...")
    finally:
        socket_io_thread.shutdown.set()
        socket_io_thread.join()


if __name__ == "__main__":
    if sys.platform == "linux":
        from ctypes import CDLL

        libc = CDLL("libc.so.6")
        PR_SET_NAME = 15  # https://github.com/torvalds/linux/blob/2df0c02dab829dd89360d98a8a1abaa026ef5798/include/uapi/linux/prctl.h#L56
        libc.prctl(PR_SET_NAME, b"agent")

    asyncio.run(main())
