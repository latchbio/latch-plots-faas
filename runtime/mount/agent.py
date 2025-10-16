import asyncio
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
from typing import Literal

import anthropic
from anthropic.types import MessageParam, ToolParam
from config_loader import build_full_instruction
from lplots import _inject
from pydantic import BaseModel
from socketio_thread import SocketIoThread
from utils import auth_token_sdk, gql_query, nucleus_url, pod_id

sandbox_root = os.environ.get("LATCH_SANDBOX_ROOT")
if sandbox_root:
    import pathlib
    original_path_new = pathlib.Path.__new__

    def patched_path_new(cls, *args, **kwargs):
        if args and args[0] == "/root/.latch":
            return original_path_new(cls, sandbox_root, *args[1:], **kwargs)
        return original_path_new(cls, *args, **kwargs)

    pathlib.Path.__new__ = patched_path_new

AGENT_DEBUG = True


class Mode(Enum):
    planning = "planning"
    executing = "executing"
    debugging = "debugging"


class PlanItem(BaseModel):
    id: str
    description: str
    status: Literal["todo", "in_progress", "done"]


class PlanDiff(BaseModel):
    action: Literal["add", "update", "complete"]
    id: str
    description: str


class NotebookResponse(BaseModel):
    plan: list[PlanItem]
    plan_diff: list[PlanDiff]
    summary: list[str] | None = None
    questions: list[str] | None = None
    next_status: Literal["executing", "fixing", "thinking", "awaiting_user_response", "awaiting_cell_execution", "done"]


from typing_extensions import TypedDict


class PlanItemPayload(TypedDict):
    id: str
    description: str
    status: Literal["todo", "in_progress", "done"]


class PlanDiffPayload(TypedDict):
    action: Literal["add", "update", "complete"]
    id: str
    description: str


"""Deprecated file-based history types removed."""


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
    instructions_context: str = ""
    current_structured_output: NotebookResponse | None = None
    current_request_id: str | None = None
    should_auto_continue: bool = False

    pending_messages: asyncio.Queue = field(default_factory=asyncio.Queue)
    conversation_task: asyncio.Task | None = None
    conversation_running: bool = False
    system_prompt: str = ""
    agent_session_id: int | None = None
    pending_action_context: dict[str, dict] = field(default_factory=dict)

    mode_config: dict[Mode, tuple[str, int | None]] = field(default_factory=lambda: {
        Mode.planning: ("claude-sonnet-4-5-20250929", 4096),
        Mode.executing: ("claude-sonnet-4-5-20250929", 1024),
        Mode.debugging: ("claude-sonnet-4-5-20250929", 2048),
    })

    async def send(self, msg: dict[str, object]) -> None:
        msg_type = msg.get("type", "unknown")
        print(f"[agent] Sending message: {msg_type}", flush=True)
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
                if role in ("user", "assistant") and (isinstance(content, (str, list))):
                    anthropic_messages.append({"role": role, "content": content})

        return anthropic_messages

    async def _insert_history(self, *, event_type: str, payload: dict, request_id: str | None = None, tx_id: str | None = None) -> None:
        if self.agent_session_id is None:
            raise RuntimeError("[insert history] Agent session ID is not set")

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
        if self.agent_session_id is None:
            raise RuntimeError("[mark all history removed] Agent session ID is not set")

        resp = await gql_query(
            auth=auth_token_sdk,
            query="""
                query AgentHistoryIds($sessionId: BigInt!) {
                    agentHistories(condition: {sessionId: $sessionId, removed: false}) {
                        nodes { id }
                    }
                }
            """,
            variables={"sessionId": str(self.agent_session_id)},
        )
        nodes = resp.get("data", {}).get("agentHistories", {}).get("nodes", [])

        # todo(aidan): move to function to batch query
        for n in nodes:
            hid = n.get("id")
            if hid is None:
                continue

            await gql_query(
                auth=auth_token_sdk,
                query="""
                    mutation RemoveHistory($id: BigInt!) {
                        updateAgentHistory(input: {id: $id, patch: {removed: true}}) {
                            clientMutationId
                        }
                    }
                """,
                variables={"id": str(hid)},
            )
        await self._notify_history_updated()

    async def atomic_operation(self, action: str, params: dict) -> dict:
        self.operation_counter += 1

        if self.mode == Mode.planning and action in {"create_cell", "edit_cell", "run_cell", "delete_cell"}:
            self.set_mode(Mode.executing)

        tx_id = f"tx_{uuid.uuid4().hex[:12]}"
        loop = asyncio.get_running_loop()
        response_future = loop.create_future()
        self.pending_operations[tx_id] = response_future

        try:
            if AGENT_DEBUG:
                print(f"[agent] -> {action}")
            # Track action+params for history once we receive the response
            self.pending_action_context[tx_id] = {"action": action, "params": params}
            await self.send({"type": "agent_action", "action": action, "params": params, "tx_id": tx_id})
        except Exception as e:
            self.pending_operations.pop(tx_id, None)
            return {"status": "error", "error": f"Send failed: {e!s}"}

        try:
            result = await asyncio.wait_for(response_future, timeout=10.0)
            return result
        except TimeoutError:
            return {"status": "error", "error": "Operation timeout", "tx_id": tx_id}
        finally:
            self.pending_operations.pop(tx_id, None)

    async def handle_action_response(self, msg: dict[str, object]) -> None:
        tx_id = msg.get("tx_id")
        fut = self.pending_operations.get(tx_id)
        if fut and not fut.done():
            fut.set_result(msg)

        # Persist a normalized agent_action history event mirroring the previous
        # frontend behavior, when possible.
        ctx = self.pending_action_context.pop(tx_id, None)
        if ctx is not None:
            try:
                action_name = ctx.get("action")
                params = ctx.get("params", {})
                status = msg.get("status")
                if status != "success":
                    return

                payload_action: dict | None = None
                if action_name == "create_cell":
                    payload_action = {
                        "task": "create_cell",
                        "cell_id": msg.get("cell_id"),
                        "tf_id": msg.get("tf_id"),
                        "auto_run": params.get("auto_run", False),
                        "source": params.get("source"),
                        "position": params.get("position"),
                        "cell_name": params.get("title"),
                    }
                elif action_name == "create_markdown_cell":
                    payload_action = {
                        "task": "create_markdown_cell",
                        "cell_id": msg.get("cell_id"),
                        "source": params.get("source"),
                    }
                elif action_name == "edit_cell":
                    payload_action = {
                        "task": "edit_cell",
                        "cell_id": params.get("cell_id"),
                        "source": params.get("source"),
                    }
                elif action_name == "run_cell":
                    payload_action = {
                        "task": "run_cell",
                        "cell_id": params.get("cell_id"),
                        **({"tf_id": msg.get("tf_id")} if msg.get("tf_id") else {}),
                    }
                elif action_name == "stop_cell":
                    payload_action = {
                        "task": "stop_cell",
                        "cell_id": params.get("cell_id"),
                        **({"tf_id": msg.get("tf_id")} if msg.get("tf_id") else {}),
                    }
                elif action_name == "set_widget":
                    value = params.get("value")
                    try:
                        parsed_value = json.loads(value) if isinstance(value, str) else value
                    except Exception:
                        parsed_value = None
                    payload_action = {
                        "task": "set_widget",
                        "widget_key": params.get("key"),
                        **({"widget_value": parsed_value} if parsed_value is not None else {}),
                    }

                if payload_action is not None:
                    await self._insert_history(
                        event_type="agent_action",
                        payload={
                            "type": "agent_action",
                            "action": payload_action,
                            "timestamp": int(time.time() * 1000),
                        },
                    )
            except Exception:
                pass

    def set_mode(self, mode: Mode) -> None:
        if mode == self.mode:
            return

        self.mode = mode
        print(f"[agent] Mode changed to {mode.value}")

    def _message_has_thinking(self, msg: dict) -> bool:
        if msg.get("role") != "assistant":
            return False
        content = msg.get("content", [])
        if not isinstance(content, list):
            return False
        for block in content:
            block_type = block.get("type") if isinstance(block, dict) else getattr(block, "type", None)
            if block_type in ("thinking", "redacted_thinking"):
                return True
        return False

    async def _wait_for_message(self) -> None:
        msg = await self.pending_messages.get()

        msg_type = msg.get("type")

        # Internal no-op used to advance the loop after tool execution/results
        if msg_type == "resume":
            if AGENT_DEBUG:
                print("[agent] Resuming turn after tool results")
            return

        if msg_type == "user_query":
            self.current_request_id = msg.get("request_id")

            # Persist exact user message for reconstruction
            await self._insert_history(
                event_type="anthropic_message",
                payload={
                    "type": "anthropic_message",
                    "role": "user",
                    "content": msg["content"],
                    "timestamp": int(time.time() * 1000),
                },
                request_id=self.current_request_id,
            )

            if AGENT_DEBUG:
                content_preview = msg["content"][:100]
                print(f"[agent] User query received: {content_preview}...")

        elif msg_type == "cell_result":
            cell_id = msg["cell_id"]
            success = msg.get("success", True)

            if success:
                result_content = f"✓ Cell {cell_id} executed successfully"
                if AGENT_DEBUG:
                    print(f"[agent] Cell {cell_id} succeeded")
            else:
                exception = msg.get("exception", "Unknown error")
                result_content = f"✗ Cell {cell_id} execution failed:\n```\n{exception}\n```"
                if AGENT_DEBUG:
                    print(f"[agent] Cell {cell_id} failed")

            await self._insert_history(
                event_type="agent_action",
                payload={
                    "type": "agent_action",
                    "action": {
                        "task": "cell_result",
                        "cell_id": cell_id,
                        "success": success,
                        **({"exception": exception} if not success else {}),
                    },
                    "timestamp": int(time.time() * 1000),
                },
            )

        elif msg_type == "stop":
            self.conversation_running = False
            if AGENT_DEBUG:
                print("[agent] Stop signal received")

    async def _send_agent_result(self) -> None:
        if not self.current_request_id:
            return

        should_continue = self.should_auto_continue
        self.should_auto_continue = False

        if self.mode == Mode.executing:
            self.set_mode(Mode.planning)

        structured = self.current_structured_output.model_dump() if self.current_structured_output is not None else None
        await self._insert_history(
            event_type="agent_result",
            payload={
            "type": "agent_result",
            "request_id": self.current_request_id,
                "responses": [],
                **({"structured_output": structured} if structured else {}),
                "timestamp": int(time.time() * 1000),
            },
            request_id=self.current_request_id,
        )
        self.current_structured_output = None

        if should_continue:
            if AGENT_DEBUG:
                print("[agent] Auto-continuing as requested by model")
            await self.pending_messages.put({
                "type": "user_query",
                "content": "Continue with the next step.",
                "request_id": self.current_request_id,
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

            if AGENT_DEBUG:
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
                if AGENT_DEBUG:
                    print(f"[tool] create_cell -> {msg}")
                return msg
            return f"Failed to create cell: {result.get('error', 'Unknown error')}"

        async def create_markdown_cell(args: dict) -> str:
            position = args["position"]
            code = args["code"]

            if position < 0:
                return "Error: Position must be non-negative"

            if AGENT_DEBUG:
                code_preview = code[:60] + "..." if len(code) > 60 else code
                print(f"[tool] create_markdown_cell pos={position} code={code_preview!r}")

            params = {
                "position": position,
                "cell_type": "markdown",
                "source": code,
            }

            result = await self.atomic_operation("create_markdown_cell", params)
            if result.get("status") == "success":
                cell_id = result.get("cell_id", "unknown")
                msg = f"Created markdown cell at position {position} (ID: {cell_id})"
                if AGENT_DEBUG:
                    print(f"[tool] create_markdown_cell -> {msg}")
                return msg
            return f"Failed to create cell: {result.get('error', 'Unknown error')}"

        async def edit_cell(args: dict) -> str:
            cell_id = args["cell_id"]
            new_code = args["new_code"]

            if AGENT_DEBUG:
                print(f"[tool] edit_cell id={cell_id}")

            params = {
                "cell_id": cell_id,
                "source": new_code,
                "auto_run": True
            }

            result = await self.atomic_operation("edit_cell", params)
            if result.get("status") == "success":
                msg = f"Cell {cell_id} edited successfully"
                if AGENT_DEBUG:
                    print(f"[tool] edit_cell -> {msg}")
                return msg
            return f"Failed to edit cell: {result.get('error', 'Unknown error')}"

        async def delete_cell(args: dict) -> str:
            cell_id = args["cell_id"]

            if AGENT_DEBUG:
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
                if AGENT_DEBUG:
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
            context_result = await self.atomic_operation("get_context", {})
            if context_result.get("status") != "success":
                return "Failed to get notebook context"

            cells = context_result.get("context", {}).get("cells", [])
            deleted_count = 0

            for cell in reversed(cells):
                cell_id = cell.get("cell_id")
                if cell_id:
                    result = await self.atomic_operation("delete_cell", {"cell_id": cell_id})
                    if result.get("status") == "success":
                        deleted_count += 1

            return f"Deleted {deleted_count} cells from the notebook"

        async def get_notebook_context(args: dict) -> str:
            params = {}

            result = await self.atomic_operation("get_context", params)
            if result.get("status") != "success":
                return f"Failed to get context: {result.get('error', 'Unknown error')}"

            context = result.get("context", {})
            cell_count = context.get("cell_count", 0)
            cells = context.get("cells", [])

            summary = f"Notebook has {cell_count} cell(s):\n"
            for cell in cells:
                index = cell.get("index", "?")
                cell_id = cell.get("cell_id", "?")
                cell_type = cell.get("cell_type", "unknown")
                status = cell.get("status", "idle")
                source = cell.get("source", "")
                tf_id = cell.get("tf_id", "?")

                source_preview = source[:500] + "..." if len(source) > 500 else source
                source_preview = source_preview.replace("\n", " ")

                summary += f"\n[{index}] ({cell_type}, {status}, cell_id: {cell_id}, tf_id: {tf_id})"
                if source_preview:
                    summary += f": {source_preview}"

                widget_summary = self._format_widget_summaries(cell.get("widgets") or [])
                if widget_summary:
                    summary += f"\n  Widgets: {widget_summary}"

            return summary

        async def set_widget(args: dict) -> str:
            key = args.get("key")
            if not key:
                return "Widget key is required"

            value = args.get("value")
            if value is None:
                return "Widget value is required"

            if AGENT_DEBUG:
                print(f"[tool] set_widget key={key} value={value!r}")

            params = {"key": key, "value": json.dumps(value)}

            result = await self.atomic_operation("set_widget", params)

            if result.get("status") == "success":
                return f"Updated widget value for: {key}"

            return f"Failed to update widget value: {result.get('error', 'Unknown error')}"

        async def submit_response(args: dict) -> str:
            try:
                summary = args.get("summary")
                if not isinstance(summary, list):
                    summary = None

                questions = args.get("questions")
                if not isinstance(questions, list):
                    questions = None

                next_status = args.get("next_status")
                if not isinstance(next_status, str) or next_status not in {"executing", "fixing", "thinking", "awaiting_user_response", "awaiting_cell_execution", "done"}:
                    if AGENT_DEBUG:
                        print(f"[agent] Invalid next_status: {next_status}")
                    return "Please provide a valid next_status"

                should_continue = args.get("continue", False)

                plan_items = args.get("plan", [])
                plan_diff_items = args.get("plan_diff", [])

                print("[tool] submit_response called with:")
                print(f"  - plan: {len(plan_items)} items")
                for item in plan_items:
                    print(f"    - [{item.get('status')}] {item.get('id')}: {item.get('description')}")
                print(f"  - plan_diff: {len(plan_diff_items)} items")
                for diff in plan_diff_items:
                    print(f"    - [{diff.get('action')}] {diff.get('id')}: {diff.get('description')}")
                print(f"  - continue: {should_continue}")

                self.current_structured_output = NotebookResponse(
                    plan=[PlanItem(**item) for item in plan_items],
                    plan_diff=[PlanDiff(**item) for item in plan_diff_items],
                    summary=summary,
                    questions=questions,
                    next_status=next_status,
                )

                if should_continue and self.executing_cells:
                    print(f"[tool] Cannot continue - {len(self.executing_cells)} cells still executing: {self.executing_cells}")
                    self.should_auto_continue = False
                else:
                    self.should_auto_continue = should_continue

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
            "name": "get_notebook_context",
            "description": "Get the current state of the notebook including all cells and their content.",
            "input_schema": {
                "type": "object",
                "properties": {},
            },
        })
        self.tool_map["get_notebook_context"] = get_notebook_context

        self.tools.append({
            "name": "submit_response",
            "description": "Submit the final response with plan, plan_diff, summary, next_status, and questions. Call this at the end of every turn.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "plan": {"type": "array", "description": "List of plan items"},
                    "plan_diff": {"type": "array", "description": "List of plan diff items"},
                    "summary": {"type": "array", "description": "List of summary bullet points or null"},
                    "questions": {"type": "array", "description": "List of questions for the user or null"},
                    "next_status": {"type": "string", "description": "What the agent will do next", "enum": ["executing", "fixing", "thinking", "awaiting_user_response", "awaiting_cell_execution", "awaiting_user_widget_input", "done"]},
                    "continue": {
                        "type": "boolean",
                        "description": "Set to true to immediately continue to the next step without waiting for user input. Set to false when waiting for user input or when all work is complete.",
                        "default": False
                    },
                },
                "required": ["plan", "plan_diff", "summary", "questions"],
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

    async def run_agent_loop(self) -> None:
        assert self.client is not None, "Client not initialized"

        self.conversation_running = True
        turn = 0

        while self.conversation_running:
            await self._wait_for_message()
            if not self.conversation_running:
                break

            turn += 1

            model, thinking_budget = self.mode_config.get(self.mode, ("claude-sonnet-4-5-20250929", 1024))

            if AGENT_DEBUG:
                print(f"[agent] Turn {turn}, mode={self.mode}, thinking_budget={thinking_budget}")

            if thinking_budget is not None:
                max_tokens = thinking_budget + 4096
            else:
                max_tokens = 4096

            def _has_content(m: dict) -> bool:
                c = m.get("content")
                if isinstance(c, str):
                    return bool(c.strip())
                if isinstance(c, list):
                    return len(c) > 0
                return False

            api_messages = await self._build_messages_from_db()

            kwargs = {
                "model": model,
                "max_tokens": max_tokens,
                "system": self.system_prompt,
                "messages": api_messages,
                "tools": self.tools,
            }
            use_beta_api = False
            if thinking_budget is not None:
                kwargs["thinking"] = {
                    "type": "enabled",
                    "budget_tokens": thinking_budget,
                }
                kwargs["betas"] = ["interleaved-thinking-2025-05-14"]
                use_beta_api = True

            try:
                start_time = time.process_time()

                if use_beta_api:
                    response = await self.client.beta.messages.create(**kwargs)
                else:
                    response = await self.client.messages.create(**kwargs)

                duration = time.process_time() - start_time

            except Exception as e:
                print(f"[agent] API error: {e}", flush=True)
                await self.send({
                    "type": "agent_error",
                    "error": f"API error: {e!s}",
                    "fatal": False
                })
                continue

            # Persist exact anthropic response so future turns can be reconstructed without
            # lossy event mapping.
            try:
                await self._insert_history(
                    event_type="anthropic_message",
                    payload={
                        "type": "anthropic_message",
                        "role": "assistant",
                        "content": response.content,
                        "timestamp": int(time.time() * 1000),
                    },
                )
            except Exception:
                # Non-fatal; continue with event-style persistence below
                pass

            for block in response.content:
                if isinstance(block, dict):
                    block_type = block.get("type")
                else:
                    block_type = getattr(block, "type", None)

                if block_type in ("thinking", "redacted_thinking"):
                    thinking_text = block.get("thinking") if isinstance(block, dict) else getattr(block, "thinking", None)
                    if thinking_text:
                        if AGENT_DEBUG:
                            print(f"[agent] Thinking:\n{thinking_text}")
                        await self._insert_history(
                            event_type="agent_thinking",
                            payload={
                            "type": "agent_thinking",
                            "thoughts": thinking_text,
                            "duration": duration,
                                "timestamp": int(time.time() * 1000),
                            },
                        )
                    else:
                        print("[agent] Thinking block present (redacted)")

            if response.stop_reason == "end_turn":
                if AGENT_DEBUG:
                    print("[agent] Turn ended without submit_response; emitting agent_result to close the turn")
                await self._send_agent_result()
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

                        if AGENT_DEBUG:
                            print(f"[agent] Executing tool: {tool_name} (id={tool_id})")

                        handler = self.tool_map.get(tool_name)
                        if handler:
                            try:
                                result = await handler(tool_input)
                                tool_results.append({
                                    "type": "tool_result",
                                    "tool_use_id": tool_id,
                                    "content": result,
                                })
                            except Exception as e:
                                print(f"[agent] Tool error: {tool_name}: {e}", flush=True)
                                tool_results.append({
                                    "type": "tool_result",
                                    "tool_use_id": tool_id,
                                    "content": f"Error executing tool: {e!s}",
                                    "is_error": True,
                                })
                        else:
                            tool_results.append({
                                "type": "tool_result",
                                "tool_use_id": tool_id,
                                "content": f"Unknown tool: {tool_name}",
                                "is_error": True,
                            })

                if tool_results:
                    for block in response.content:
                        block_type = block.get("type") if isinstance(block, dict) else getattr(block, "type", None)
                        if block_type == "tool_use":
                            tool_id = block.get("id") if isinstance(block, dict) else block.id
                            tool_name = block.get("name") if isinstance(block, dict) else block.name
                            tool_input = block.get("input") if isinstance(block, dict) else block.input
                            await self._insert_history(
                                event_type="agent_action",
                                payload={
                                    "type": "agent_action",
                                    "action": {
                                        "task": "tool_use",
                                        "id": tool_id,
                                        "name": tool_name,
                                        "input": tool_input,
                                    },
                                    "timestamp": int(time.time() * 1000),
                                },
                            )
                    # Persist the exact user tool_result message shape too
                    try:
                        await self._insert_history(
                            event_type="anthropic_message",
                            payload={
                                "type": "anthropic_message",
                                "role": "user",
                                "content": tool_results,
                                "timestamp": int(time.time() * 1000),
                            },
                        )
                    except Exception:
                        pass
                    for tr in tool_results:
                        await self._insert_history(
                            event_type="agent_action",
                            payload={
                                "type": "agent_action",
                                "action": {
                                    "task": "tool_result",
                                    "tool_use_id": tr.get("tool_use_id"),
                                    "content": tr.get("content"),
                                    "is_error": tr.get("is_error", False),
                                },
                                "timestamp": int(time.time() * 1000),
                            },
                        )
                elif AGENT_DEBUG:
                    print("[agent] No tool results")

                if called_submit_response:
                    if AGENT_DEBUG:
                        print("[agent] submit_response called, sending agent_result")
                    await self._send_agent_result()
                else:
                    # Continue the turn by prompting the loop to run another
                    # round with the updated DB-backed history (tool use + results)
                    await self.pending_messages.put({"type": "resume"})
            elif response.stop_reason == "max_tokens":
                print("[agent] Hit max tokens", flush=True)
                await self._send_agent_result()
            else:
                if AGENT_DEBUG:
                    print(f"[agent] Unknown stop reason: {response.stop_reason}")
                await self._send_agent_result()

    async def handle_init(self, msg: dict[str, object]) -> None:
        print("[agent] Initializing", flush=True)

        try:
            context = msg.get("context", "")
            self.instructions_context = context
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

            self.system_prompt = build_full_instruction(self.instructions_context)

            self.initialized = True
            await self.send({
                "type": "agent_status",
                "status": "ready"
            })
            print("[agent] Initialization complete", flush=True)

            self.conversation_task = asyncio.create_task(self.run_agent_loop())
        except Exception as e:
            await self.send({
                "type": "agent_error",
                "error": f"Failed to initialize: {e!s}",
                "fatal": True
            })

    async def handle_query(self, msg: dict[str, object]) -> None:
        query = msg.get("query", "")
        request_id = msg.get("request_id")

        print(f"[agent] Processing query: {query}...")

        await self.pending_messages.put({
            "type": "user_query",
            "content": query,
            "request_id": request_id,
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
        print("[agent] Clearing conversation history")
        await self._mark_all_history_removed()

    async def accept(self) -> None:
        msg = await self.conn.recv()
        msg_type = msg.get("type")

        print(f"[agent] Received message: {msg_type}", flush=True)

        if msg_type == "init":
            if AGENT_DEBUG:
                print(f"[agent] Message: {msg_type}")
            await self.handle_init(msg)
        elif msg_type == "agent_query":
            if AGENT_DEBUG:
                query = msg.get("query", "")
                query_preview = query[:60] + "..." if len(query) > 60 else query
                print(f"[agent] Query: {query_preview}")
            await self.handle_query(msg)
        elif msg_type == "agent_cancel":
            if AGENT_DEBUG:
                request_id = msg.get("request_id", "unknown")
                print(f"[agent] Cancel: {request_id}")
            await self.handle_cancel(msg)
        elif msg_type == "agent_clear_history":
            if AGENT_DEBUG:
                print("[agent] Clear history request")
            await self.handle_clear_history()
        elif msg_type == "agent_action_response":
            if AGENT_DEBUG:
                action = msg.get("action", "unknown")
                status = msg.get("status", "unknown")
                print(f"[agent] {action} -> {status}")
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
    print(f"{datetime.now().isoformat()} [agent] Starting", flush=True)

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
