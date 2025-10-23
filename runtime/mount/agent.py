import asyncio
import contextlib
import json
import os
import socket
import sys
import time
import traceback
import uuid
from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum
from textwrap import dedent

from agent_utils.auto_install import anthropic
from anthropic.types import MessageParam, ToolParam
from anthropic.types.beta.beta_message import BetaMessage
from anthropic.types.message import Message
from config_loader import build_system_prompt
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
    pending_operations: dict[str, asyncio.Future] = field(default_factory=dict)
    executing_cells: set[str] = field(default_factory=set)
    tools: list[ToolParam] = field(default_factory=list)
    tool_map: dict[str, Callable] = field(default_factory=dict)
    operation_counter: int = 0
    current_request_id: str | None = None
    should_auto_continue: bool = False
    pending_auto_continue: bool = False

    pending_messages: asyncio.Queue = field(default_factory=asyncio.Queue)
    conversation_task: asyncio.Task | None = None
    conversation_running: bool = False
    system_prompt: str = ""
    agent_session_id: int | None = None

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
                        nodes { id payload }
                    }
                }
            """,
            variables={"sessionId": str(self.agent_session_id)},
        )
        nodes = resp.get("data", {}).get("agentHistories", {}).get("nodes", [])
        return [n.get("payload") for n in nodes if n.get("payload")]

    async def _build_messages_from_db(self) -> list[MessageParam]:
        history = await self._fetch_history_from_db()
        anthropic_messages: list[MessageParam] = []

        for item in history:
            t = item.get("type") if isinstance(item, dict) else None
            if t == "anthropic_message":
                role = item.get("role")
                content = item.get("content")
                if role in {"user", "assistant"} and (isinstance(content, (str, list))):
                    anthropic_messages.append({"role": role, "content": content})

        print(f"[agent] Built {len(anthropic_messages)} messages from DB")

        return anthropic_messages

    async def _insert_history(self, *, event_type: str, payload: dict, request_id: str | None = None, tx_id: str | None = None) -> None:
        assert self.agent_session_id is not None

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
        msg = await self.pending_messages.get()

        msg_type = msg.get("type")
        print(f"[agent] _wait_for_message: got message type={msg_type}")

        if msg_type == "resume":
            print("[agent] Resuming turn after tool results")
            return

        if msg_type == "user_query":
            self.current_request_id = msg.get("request_id")

            payload = {
                "type": "anthropic_message",
                "role": "user",
                "content": msg["content"],
                "timestamp": int(time.time() * 1000),
            }

            if msg.get("display_query") is not None:
                payload["display_query"] = msg["display_query"]
            if msg.get("display_nodes") is not None:
                payload["display_nodes"] = msg["display_nodes"]
            if msg.get("hidden") is not None:
                payload["hidden"] = msg["hidden"]

            await self._insert_history(
                event_type="anthropic_message",
                payload=payload,
                request_id=self.current_request_id,
            )

        elif msg_type == "cell_result":
            cell_id = msg["cell_id"]
            success = msg.get("success", True)

            if success:
                result_content = f"✓ Cell {cell_id} executed successfully"
                print(f"[agent] Cell {cell_id} succeeded")
            else:
                exception = msg.get("exception", "Unknown error")
                result_content = f"✗ Cell {cell_id} execution failed:\n```\n{exception}\n```"
                print(f"[agent] Cell {cell_id} failed")

            await self._insert_history(
                event_type="anthropic_message",
                payload={
                    "type": "anthropic_message",
                    "role": "user",
                    "content": result_content,
                    "timestamp": int(time.time() * 1000),
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

        elif msg_type == "stop":
            self.conversation_running = False
            print("[agent] Stop signal received")

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

        async def create_cell(args: dict) -> str:
            position = args["position"]
            code = args["code"]
            title = args["title"]

            if position < 0:
                return "Error: Position must be non-negative"

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
                return msg
            return f"Failed to create cell: {result.get('error', 'Unknown error')}"

        async def create_markdown_cell(args: dict) -> str:
            position = args["position"]
            code = args["code"]

            if position < 0:
                return "Error: Position must be non-negative"

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
                return msg
            return f"Failed to create cell: {result.get('error', 'Unknown error')}"

        async def edit_cell(args: dict) -> str:
            cell_id = args["cell_id"]
            new_code = args["new_code"]

            print(f"[tool] edit_cell id={cell_id}")

            params = {
                "cell_id": cell_id,
                "source": new_code,
                "auto_run": True
            }

            result = await self.atomic_operation("edit_cell", params)
            if result.get("status") == "success":
                msg = f"Cell {cell_id} edited successfully"
                print(f"[tool] edit_cell -> {msg}")
                return msg
            return f"Failed to edit cell: {result.get('error', 'Unknown error')}"

        async def delete_cell(args: dict) -> str:
            cell_id = args["cell_id"]

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
                return msg
            return f"Failed to delete cell: {result.get('error', 'Unknown error')}"

        async def run_cell(args: dict) -> str:
            cell_id = args["cell_id"]
            params = {"cell_id": cell_id}

            await self.send({
                "type": "agent_action",
                "action": "run_cell",
                "params": params,
            })
            self.executing_cells.add(cell_id)

            return f"Cell {cell_id} execution started"

        async def stop_cell(args: dict) -> str:
            cell_id = args["cell_id"]
            params = {"cell_id": cell_id}

            result = await self.atomic_operation("stop_cell", params)
            if result.get("status") == "success":
                self.executing_cells.discard(cell_id)
                return f"Stopped cell {cell_id}"
            return f"Failed to stop cell {cell_id}: {result.get('error', 'Unknown error')}"

        async def delete_all_cells(args: dict) -> str:
            context_result = await self.atomic_operation("get_context")
            if context_result.get("status") != "success":
                error_msg = context_result.get("error", "Unknown error")
                return f"Failed to delete cells: {error_msg}"

            cells = context_result.get("context", {}).get("cells", [])
            deleted_count = 0

            for cell in reversed(cells):
                cell_id = cell.get("cell_id")
                if cell_id:
                    result = await self.atomic_operation("delete_cell", {"cell_id": cell_id})
                    if result.get("status") == "success":
                        deleted_count += 1

            return f"Deleted {deleted_count} cells from the notebook"

        async def set_widget(args: dict) -> str:
            key = args.get("key")
            if not key:
                return "Widget key is required"

            value = args.get("value")
            if value is None:
                return "Widget value is required"

            print(f"[tool] set_widget key={key} value={value!r}")

            params = {"key": key, "value": json.dumps(value)}

            result = await self.atomic_operation("set_widget", params)

            if result.get("status") == "success":
                return f"Updated widget value for: {key}"

            return f"Failed to update widget value: {result.get('error', 'Unknown error')}"

        async def submit_response(args: dict) -> str:
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
                    return "Please provide a valid next_status"

                should_continue = args.get("continue", False)

                plan_items = args.get("plan", [])
                plan_diff_items = args.get("plan_diff", [])

                print("[tool] submit_response called with:")
                print(f"  - next_status: {next_status}")
                print(f"  - plan: {len(plan_items)} items")
                for item in plan_items:
                    print(f"    - [{item.get('status')}] {item.get('id')}: {item.get('description')}")
                print(f"  - plan_diff: {len(plan_diff_items)} items")
                for diff in plan_diff_items:
                    print(f"    - [{diff.get('action')}] {diff.get('id')}: {diff.get('description')}")
                print(f"  - summary: {summary}")
                print(f"  - questions: {questions}")
                print(f"  - continue: {should_continue}")

                if should_continue and self.executing_cells:
                    print(f"[tool] Deferring auto-continue - {len(self.executing_cells)} cells still executing: {self.executing_cells}")
                    self.should_auto_continue = False
                    self.pending_auto_continue = True
                else:
                    self.should_auto_continue = should_continue
                    self.pending_auto_continue = False

                return "Response submitted successfully"
            except Exception as e:
                print(f"[tool] submit_response error: {e}")
                import traceback
                traceback.print_exc()
                return f"Error submitting response: {e!s}"

        self.tools.append({
            "name": "create_cell",
            "description": "Create a new code cell at specified position. The cell will automatically run after creation.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "position": {"type": "integer", "description": "Position to insert the cell"},
                    "code": {"type": "string", "description": "Python code for the cell"},
                    "title": {"type": "string", "description": "Title for the cell"},
                },
                "required": ["position", "code", "title"],
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
                },
                "required": ["position", "code"],
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
                },
                "required": ["cell_id", "new_code"],
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
                },
                "required": ["cell_id"],
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
                },
                "required": ["cell_id"],
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
                },
                "required": ["cell_id"],
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
            "name": "submit_response",
            "description": "Submit the final response with plan, plan_diff, next_status, questions, and an optional summary. Call this at the end of every turn.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "plan": {"type": "array", "description": "List of plan items"},
                    "plan_diff": {"type": "array", "description": "List of plan diff items"},
                    "summary": {"type": "string", "description": "Optional summary text describing what was accomplished. Use markdown formatting with bullet points if needed. Omit if no summary needed."},
                    "questions": {"type": "string", "description": "Optional question text for the user. Omit if no questions needed."},
                    "next_status": {"type": "string", "description": "What the agent will do next", "enum": ["executing", "fixing", "thinking", "awaiting_user_response", "awaiting_cell_execution", "awaiting_user_widget_input", "done"]},
                    "continue": {
                        "type": "boolean",
                        "description": "Set to true to immediately continue to the next step without waiting for user input. Set to false when waiting for user input or when all work is complete.",
                        "default": False
                    },
                },
                "required": ["plan", "plan_diff", "next_status"],
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
                        },
                        "required": ["key", "value"],
                    },
                })
        self.tool_map["set_widget"] = set_widget

    async def _get_notebook_context(self) -> str:
        context_result, globals_result, reactivity_result = await asyncio.gather(
            self.atomic_operation("get_context"),
            self.atomic_operation("request_globals_summary"),
            self.atomic_operation("request_reactivity_summary")
        )

        if context_result.get("status") != "success":
            return f"Failed to get context: {context_result.get('error', 'Unknown error')}"

        context = context_result.get("context", {})
        cell_count = context.get("cell_count", 0)
        cells = context.get("cells", [])

        globals_data: dict[str, object] | None = None
        if globals_result.get("status") == "success":
            globals_data = globals_result.get("summary", {})

        reactivity_summary: str | None = None
        if reactivity_result.get("status") == "success":
            reactivity_summary = reactivity_result.get("summary")

        sections = []

        if reactivity_summary is not None:
            sections.append(reactivity_summary)

        cell_section = ["## Notebook Cells", f"\nTotal: {cell_count} cells\n"]

        for cell in cells:
            index = cell.get("index", "?")
            cell_id = cell.get("cell_id", "?")
            cell_type = cell.get("cell_type", "unknown")
            source = cell.get("source", "")
            status = cell.get("status", "idle")
            tf_id = cell.get("tf_id", "?")

            cell_section.extend((f"\n### Cell [{index}] (ID: {cell_id}, TF ID: {tf_id})", f"Type: {cell_type}", f"Status: {status}"))

            if source is not None:
                source_display = source[:800] if len(source) > 800 else source
                truncated = " [TRUNCATED]" if len(source) > 800 else ""
                cell_section.append(f"```python\n{source_display}\n```{truncated}")

            widgets = cell.get("widgets", None)
            if widgets is not None:
                widget_strs = []
                for w in widgets:
                    w_type = w.get("type", "unknown")
                    w_key = w.get("key", "")
                    w_label = w.get("label", "")
                    if w_label:
                        widget_strs.append(f"{w_type} ({w_label}) [{w_key}]")
                    else:
                        widget_strs.append(f"{w_type} [{w_key}]")
                if widget_strs:
                    cell_section.append(f"Widgets: {', '.join(widget_strs)}")

        sections.append("\n".join(cell_section))

        if globals_data is not None and len(globals_data) > 0:
            global_section = ["\n## Global Variables", f"\nTotal: {len(globals_data)} variables\n"]

            for var_name in sorted(globals_data.keys()):
                var_info = globals_data[var_name]
                if isinstance(var_info, dict):
                    var_type = var_info.get("type", "unknown")
                    global_section.append(f"\n**{var_name}** ({var_type})")
                    for key, value in var_info.items():
                        if key != "type":
                            global_section.append(f"  - {key}: {value}")
                else:
                    global_section.append(f"\n**{var_name}**: {var_info}")

            sections.append("\n".join(global_section))

        return "\n\n".join(sections)

    async def run_agent_loop(self) -> None:
        assert self.client is not None, "Client not initialized"

        self.conversation_running = True
        turn = 0

        while self.conversation_running:
            await self._wait_for_message()
            if not self.conversation_running:
                break

            api_messages = await self._build_messages_from_db()
            if not api_messages or api_messages[-1].get("role") != "user":
                continue

            turn += 1

            model, thinking_budget = self.mode_config.get(self.mode, ("claude-sonnet-4-5-20250929", 1024))

            print(f"[agent] Turn {turn}, mode={self.mode}, thinking_budget={thinking_budget}")

            notebook_context = await self._get_notebook_context()

            print(f"[agent] Notebook context:\n{notebook_context}")

            system_prompt_with_context = dedent(
                f"""
                {self.system_prompt}
                <notebook_context>
                {notebook_context}
                </notebook_context>
                """
            )

            if thinking_budget is not None:
                max_tokens = thinking_budget + 4096
            else:
                max_tokens = 4096

            can_use_thinking = True
            if thinking_budget is not None and len(api_messages) > 0:
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

            kwargs = {
                "model": model,
                "max_tokens": max_tokens,
                "system": system_prompt_with_context,
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
                start_time = time.time()

                if use_beta_api:
                    response: BetaMessage = await self.client.beta.messages.create(**kwargs)
                else:
                    response: Message = await self.client.messages.create(**kwargs)

                duration_seconds = time.time() - start_time

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
                    event_type="anthropic_message",
                    payload={
                        "type": "anthropic_message",
                        "role": "assistant",
                        "content": response_content,
                        "timestamp": int(time.time() * 1000),
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
                        if handler:
                            try:
                                result = await handler(tool_input)
                                tool_results.append({
                                    "type": "tool_result",
                                    "tool_use_id": tool_id,
                                    "content": json.dumps({
                                        "summary": result,
                                        "error": None,
                                    }),
                                })
                            except Exception as e:
                                print(f"[agent] Tool error: {tool_name}: {e}")
                                tool_results.append({
                                    "type": "tool_result",
                                    "tool_use_id": tool_id,
                                    "content": json.dumps({
                                        "summary": None,
                                        "error": f"Error executing tool: {e!s}",
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
                    await self._insert_history(
                        event_type="anthropic_message",
                        payload={
                            "type": "anthropic_message",
                            "role": "user",
                            "content": tool_results,
                            "timestamp": int(time.time() * 1000),
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

            self.system_prompt = build_system_prompt()

            messages = await self._build_messages_from_db()
            tool_use_ids = set()
            tool_result_ids = set()

            print(f"[agent] Checking for pending tools in {len(messages)} messages")

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
                            print(f"[agent]   Found tool_use: id={tool_id}, name={block.get('name')}")
                    elif block.get("type") == "tool_result":
                        tool_use_id = block.get("tool_use_id")
                        if tool_use_id:
                            tool_result_ids.add(tool_use_id)
                            print(f"[agent]   Found tool_result: tool_use_id={tool_use_id}")

            pending_tool_ids = tool_use_ids - tool_result_ids

            print(f"[agent] Tool use count: {len(tool_use_ids)}, Tool result count: {len(tool_result_ids)}, Pending: {len(pending_tool_ids)}")

            if pending_tool_ids:
                print(f"[agent] Found {len(pending_tool_ids)} pending tool calls, closing them")

                await self._insert_history(
                    event_type="anthropic_message",
                    payload={
                        "type": "anthropic_message",
                        "role": "user",
                        "content": [
                            {
                                "type": "tool_result",
                                "tool_use_id": tool_id,
                                "content": json.dumps({
                                    "summary": None,
                                    "error": "Tool call cancelled - session was ended before completion",
                                }),
                            }
                            for tool_id in pending_tool_ids
                        ],
                        "timestamp": int(time.time() * 1000),
                    },
                )
                await self._notify_history_updated()

            self.initialized = True
            await self.send({
                "type": "agent_status",
                "status": "ready"
            })
            print("[agent] Initialization complete")

            self._start_conversation_loop()

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
                    print(f"[agent] Reconnected while {next_status}, prompting LLM to check state and retry")
                    await self.pending_messages.put({
                        "type": "user_query",
                        "content": "The session was reconnected. You were waiting for an action to complete, but it may have finished or failed while offline. Please retry any incomplete actions or continue with your plan.",
                        "request_id": None,
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

        print(f"[agent] Processing query: {query}...")

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

    async def handle_cancel(self, msg: dict[str, object]) -> None:
        request_id = msg.get("request_id", "unknown")
        print(f"[agent] Cancelling request {request_id}")

        self.conversation_running = False

        if self.conversation_task and not self.conversation_task.done():
            self.conversation_task.cancel()
            try:
                await self.conversation_task
            except asyncio.CancelledError:
                pass

    def _format_widget_summaries(self, widgets: list[dict]) -> str:
        if not widgets:
            return ""

        widget_summaries = []
        for widget in widgets:
            widget_key = widget.get("key") or "?"
            widget_type = widget.get("type") or "unknown"
            widget_value = widget.get("value") or ""
            widget_label = widget.get("label") or ""

            widget_desc = f"key={widget_key}, type={widget_type}"
            if widget_value:
                widget_desc += f", value={widget_value}"
            if widget_label:
                widget_desc += f", label={widget_label}"
            widget_summaries.append(widget_desc)

        return ", ".join(widget_summaries)

    async def handle_clear_history(self) -> None:
        await self._clear_running_state()
        await self._mark_all_history_removed()
        self._start_conversation_loop()

    async def accept(self) -> None:
        msg = await self.conn.recv()
        msg_type = msg.get("type")

        print(f"[agent] Received message: {msg_type}")

        if msg_type == "init":
            print(f"[agent] Message: {msg_type}")
            await self.handle_init(msg)
        elif msg_type == "agent_query":
            query = msg.get("query", "")
            query_preview = query[:60] + "..." if len(query) > 60 else query
            print(f"[agent] Query: {query_preview}")
            await self.handle_query(msg)
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
            nested_msg = msg.get("message", {})
            nested_type = nested_msg.get("type")
            if nested_type == "cell_result":
                cell_id = nested_msg.get("cell_id")
                has_exception = nested_msg.get("has_exception", False)
                exception = nested_msg.get("exception", "")

                if cell_id is not None:
                    self.executing_cells.discard(str(cell_id))

                await self.pending_messages.put({
                    "type": "cell_result",
                    "cell_id": cell_id,
                    "success": not has_exception,
                    "exception": exception,
                })
            elif nested_type == "start_cell":
                cell_id = nested_msg.get("cell_id")
                if cell_id is not None:
                    self.executing_cells.add(str(cell_id))
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
