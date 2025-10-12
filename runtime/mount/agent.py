import asyncio
import json
import os
import socket
import sys
import traceback
import uuid
from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum
from textwrap import dedent
from typing import Literal

import anthropic
from anthropic.types import MessageParam, ToolParam, ToolUseBlock
from config_loader import build_full_instruction
from lplots import _inject
from pydantic import BaseModel
from socketio_thread import SocketIoThread

sandbox_root = os.environ.get("LATCH_SANDBOX_ROOT")
if sandbox_root:
    import pathlib
    original_path_new = pathlib.Path.__new__

    def patched_path_new(cls, *args, **kwargs):
        if args and args[0] == "/root/.latch":
            return original_path_new(cls, sandbox_root, *args[1:], **kwargs)
        return original_path_new(cls, *args, **kwargs)

    pathlib.Path.__new__ = patched_path_new

AGENT_DEBUG = os.environ.get("AGENT_DEBUG") == "1"


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


from typing_extensions import TypedDict


class PlanItemPayload(TypedDict):
    id: str
    description: str
    status: Literal["todo", "in_progress", "done"]


class PlanDiffPayload(TypedDict):
    action: Literal["add", "update", "complete"]
    id: str
    description: str

class WidgetUpdatePayload(TypedDict):
    key: str
    value: str

@dataclass
class AgentHarness:
    conn: SocketIoThread
    initialized: bool = False
    api_key: str | None = None
    client: anthropic.AsyncAnthropic | None = None
    conversation_history: list[MessageParam] = field(default_factory=list)
    mode: Mode = Mode.planning
    pending_operations: dict[str, asyncio.Future] = field(default_factory=dict)
    executing_cells: set[str] = field(default_factory=set)
    tools: list[ToolParam] = field(default_factory=list)
    tool_map: dict[str, Callable] = field(default_factory=dict)
    operation_counter: int = 0
    instructions_context: str = ""
    current_structured_output: NotebookResponse | None = None
    current_request_id: str | None = None

    pending_messages: asyncio.Queue = field(default_factory=asyncio.Queue)
    conversation_task: asyncio.Task | None = None
    conversation_running: bool = False
    system_prompt: str = ""

    mode_config: dict[Mode, tuple[str, int | None]] = field(default_factory=lambda: {
        Mode.planning: ("claude-sonnet-4-5-20250929", 4096),
        Mode.executing: ("claude-sonnet-4-5-20250929", 1024),
        Mode.debugging: ("claude-sonnet-4-5-20250929", 2048),
    })

    async def send(self, msg: dict[str, object]) -> None:
        msg_type = msg.get("type", "unknown")
        print(f"[agent] Sending message: {msg_type}", flush=True)
        await self.conn.send(msg)

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

    def set_mode(self, mode: Mode) -> None:
        if mode == self.mode:
            return

        self.mode = mode
        print(f"[agent] Mode changed to {mode.value}")

    def _last_assistant_msg_has_thinking(self) -> bool:
        for msg in self.conversation_history[::-1]:
            if msg.get("role") != "assistant":
                continue

            content = msg.get("content", [])
            if not isinstance(content, list):
                return False

            has_thinking = False
            has_text = False

            for block in content:
                if isinstance(block, dict):
                    block_type = block.get("type")
                else:
                    block_type = getattr(block, "type", None)

                if block_type in ("thinking", "redacted_thinking"):
                    has_thinking = True
                elif block_type == "text":
                    has_text = True

            if has_thinking:
                return True
            if has_text:
                return False

        return False

    async def _wait_for_message(self) -> None:
        msg = await self.pending_messages.get()

        msg_type = msg.get("type")

        if msg_type == "user_query":
            self.current_request_id = msg.get("request_id")

            self.conversation_history.append({
                "role": "user",
                "content": msg["content"],
            })

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

            self.conversation_history.append({
                "role": "user",
                "content": result_content,
            })

        elif msg_type == "stop":
            self.conversation_running = False
            if AGENT_DEBUG:
                print("[agent] Stop signal received")

    async def _send_agent_result(self) -> None:
        if not self.current_request_id:
            return

        if self.mode == Mode.executing:
            self.set_mode(Mode.planning)

        text_responses = []
        for msg in reversed(self.conversation_history):
            if msg.get("role") == "user":
                break
            if msg.get("role") == "assistant":
                for block in msg.get("content", []):
                    block_type = (block.get("type") if isinstance(block, dict)
                                 else getattr(block, "type", None))
                    block_text = (block.get("text") if isinstance(block, dict)
                                 else getattr(block, "text", None))
                    if block_type == "text" and block_text:
                        text_responses.insert(0, block_text)

        response_msg = {
            "type": "agent_result",
            "status": "success",
            "responses": text_responses,
            "mode": self.mode.value,
            "request_id": self.current_request_id,
        }

        if self.current_structured_output is not None:
            response_msg["structured_output"] = self.current_structured_output.model_dump()
            self.current_structured_output = None

        await self.send(response_msg)

        self.current_request_id = None

    def init_tools(self) -> None:
        self.tools = []
        self.tool_map = {}

        async def create_cell(args: dict) -> str:
            position = args["position"]
            code = args["code"]
            title = args["title"]
            auto_run = args.get("auto_run", True)

            if position < 0:
                return "Error: Position must be non-negative"

            if AGENT_DEBUG:
                print(f'[tool] create_cell pos={position} title="{title}"')

            params = {
                "position": position,
                "cell_type": "code",
                "source": code,
                "title": title,
                "auto_run": auto_run,
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
            auto_run = args.get("auto_run", True)

            if AGENT_DEBUG:
                print(f"[tool] edit_cell id={cell_id}")

            params = {
                "cell_id": cell_id,
                "source": new_code,
                "auto_run": auto_run
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
            updates = args.get("updates")
            if not updates:
                return "No widget updates provided"

            updates_dict: dict[str, object] = {}
            payload_updates: list[WidgetUpdatePayload] = []

            if AGENT_DEBUG:
                print(f"[tool] set_widget updates={updates}")

            for item in updates:
                key = item.get("key")
                if not key:
                    return "Widget updates must include a 'key'"

                value_json = item.get("value")
                if not value_json:
                    return "Widget updates must include a 'value'"

                try:
                    parsed_value = json.loads(value_json)
                except json.JSONDecodeError:
                    parsed_value = value_json

                updates_dict[key] = parsed_value
                payload_updates.append({"key": key, "value": json.dumps(parsed_value)})

            params = {"widget_updates": payload_updates}
            result = await self.atomic_operation("set_widget", params)

            if result.get("status") == "success":
                applied_keys = ", ".join(updates_dict.keys())
                return f"Updated widget values for: {applied_keys}"

            return f"Failed to update widget values: {result.get('error', 'Unknown error')}"

        async def send_plan_update(args: dict) -> str:
            plan = args["plan"]
            plan_diff = args.get("plan_diff")

            if AGENT_DEBUG:
                print(f"[tool] send_plan_update plan_items={len(plan)} diff_items={len(plan_diff or [])}")

            try:
                plan_items = [PlanItem(**item) for item in plan]
                plan_diff_items = [PlanDiff(**item) for item in (plan_diff or [])]
            except Exception as e:
                return f"Invalid plan data: {e}"

            plan_payload = {
                "plan": [item.model_dump() for item in plan_items],
                "plan_diff": [item.model_dump() for item in plan_diff_items],
            }

            result = await self.atomic_operation("plan_update", plan_payload)
            if result.get("status") == "success":
                msg = "Plan update delivered"
                if AGENT_DEBUG:
                    print(f"[tool] send_plan_update -> {msg}")
                return msg
            return f"Failed to deliver plan update: {result.get('error', 'Unknown error')}"

        async def start_new_plan(args: dict) -> str:
            self.set_mode(Mode.planning)
            return "Started new planning session"

        async def submit_response(args: dict) -> str:
            try:
                summary = args.get("summary")
                if not isinstance(summary, list):
                    summary = None

                questions = args.get("questions")
                if not isinstance(questions, list):
                    questions = None

                self.current_structured_output = NotebookResponse(
                    plan=[PlanItem(**item) for item in args.get("plan", [])],
                    plan_diff=[PlanDiff(**item) for item in args.get("plan_diff", [])],
                    summary=summary,
                    questions=questions
                )

                if AGENT_DEBUG:
                    print(f"[tool] submit_response stored structured output")

                return "Response submitted successfully"
            except Exception as e:
                print(f"[tool] submit_response error: {e}")
                return f"Error submitting response: {e!s}"

        self.tools.append({
            "name": "create_cell",
            "description": "Create a new code cell at specified position.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "position": {"type": "integer", "description": "Position to insert the cell"},
                    "code": {"type": "string", "description": "Python code for the cell"},
                    "title": {"type": "string", "description": "Title for the cell"},
                    "auto_run": {"type": "boolean", "description": "Whether to run the cell after creation", "default": True},
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
            "description": "Replace the contents of an existing cell.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "cell_id": {"type": "string", "description": "ID of the cell to edit"},
                    "new_code": {"type": "string", "description": "New code/content for the cell"},
                    "auto_run": {"type": "boolean", "description": "Whether to run the cell after editing", "default": True},
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
            "name": "send_plan_update",
            "description": "Send plan state update to the frontend.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "plan": {"type": "array", "description": "List of plan items"},
                    "plan_diff": {"type": "array", "description": "List of plan diff items"},
                },
                "required": ["plan", "plan_diff"],
            },
        })
        self.tool_map["send_plan_update"] = send_plan_update

        self.tools.append({
            "name": "start_new_plan",
            "description": "Start a new planning session.",
            "input_schema": {
                "type": "object",
                "properties": {},
            },
        })
        self.tool_map["start_new_plan"] = start_new_plan

        self.tools.append({
            "name": "submit_response",
            "description": "Submit the final response with plan, plan_diff, summary, and questions. Call this at the end of every response.",
            "input_schema": {
                "type": "object",
                "properties": {
                    "plan": {"type": "array", "description": "List of plan items"},
                    "plan_diff": {"type": "array", "description": "List of plan diff items"},
                    "summary": {"type": "array", "description": "List of summary bullet points or null"},
                    "questions": {"type": "array", "description": "List of questions for the user or null"},
                },
                "required": ["plan", "plan_diff", "summary", "questions"],
            },
        })
        self.tool_map["submit_response"] = submit_response

        self.tools.append({
                    "name": "set_widget",
                    "description": "Set widget values by widget key.",
                    "input_schema": {
                        "type": "object",
                        "properties": {
                            "updates": {
                                "type": "array",
                                "description": "List of widget updates with keys and JSON-serializable values.",
                                "items": {
                                    "type": "object",
                                    "properties": {
                                        "key": {"type": "string", "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>"},
                                        "value": {"type": "string", "description": "JSON-serializable value string"},
                                    },
                                    "required": ["key", "value"],
                                },
                            },
                        },
                        "required": ["updates"],
                    },
                })
        self.tool_map["set_widget"] = set_widget

    async def run_agent_loop(self) -> None:
        assert self.client is not None, "Client not initialized"

        self.conversation_running = True
        turn = 0

        while self.conversation_running:
            if not self.conversation_history or self.conversation_history[-1]["role"] == "assistant":
                await self._wait_for_message()
                continue

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

            kwargs = {
                "model": model,
                "max_tokens": max_tokens,
                "system": self.system_prompt,
                "messages": self.conversation_history,
                "tools": self.tools,
            }

            use_beta_api = False
            if thinking_budget is not None and (len(self.conversation_history) == 1 or self._last_assistant_msg_has_thinking()):
                kwargs["thinking"] = {
                    "type": "enabled",
                    "budget_tokens": thinking_budget,
                }
                kwargs["betas"] = ["interleaved-thinking-2025-05-14"]
                use_beta_api = True

            try:
                if use_beta_api:
                    response = await self.client.beta.messages.create(**kwargs)
                else:
                    response = await self.client.messages.create(**kwargs)
            except Exception as e:
                print(f"[agent] API error: {e}", flush=True)
                await self.send({
                    "type": "agent_error",
                    "error": f"API error: {e!s}",
                    "fatal": False
                })
                continue

            assistant_msg = {
                "role": "assistant",
                "content": response.content,
            }
            self.conversation_history.append(assistant_msg)

            for block in response.content:
                if isinstance(block, dict):
                    block_type = block.get("type")
                else:
                    block_type = getattr(block, "type", None)

                if block_type in ("thinking", "redacted_thinking"):
                    if AGENT_DEBUG:
                        thinking_text = block.get("thinking") if isinstance(block, dict) else getattr(block, "thinking", None)
                        if thinking_text:
                            print(f"[agent] Thinking:\n{thinking_text}")
                        else:
                            print(f"[agent] Thinking block present (redacted)")

            if response.stop_reason == "end_turn":
                if AGENT_DEBUG:
                    print("[agent] Turn ended, sending agent_result")
                await self._send_agent_result()
            elif response.stop_reason == "tool_use":
                tool_results = []

                for block in response.content:
                    if isinstance(block, dict):
                        block_type = block.get("type")
                    else:
                        block_type = getattr(block, "type", None)

                    if block_type == "tool_use":
                        tool_id = block.get("id") if isinstance(block, dict) else block.id
                        tool_name = block.get("name") if isinstance(block, dict) else block.name
                        tool_input = block.get("input") if isinstance(block, dict) else block.input

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
                    self.conversation_history.append({
                        "role": "user",
                        "content": tool_results,
                    })
                else:
                    if AGENT_DEBUG:
                        print("[agent] No tool results")
            elif response.stop_reason == "max_tokens":
                print("[agent] Hit max tokens", flush=True)
                await self._send_agent_result()
            else:
                if AGENT_DEBUG:
                    print(f"[agent] Unknown stop reason: {response.stop_reason}")
                await self._send_agent_result()

    async def handle_init(self, msg: dict[str, object]) -> None:
        print("[agent] Initializing", flush=True)

        self.api_key = os.environ.get("ANTHROPIC_API_KEY")

        if not self.api_key:
            await self.send({
                "type": "agent_error",
                "error": "ANTHROPIC_API_KEY not set",
                "fatal": True
            })
            return

        try:
            context = msg.get("context", "")
            self.instructions_context = context

            self.init_tools()

            self.client = anthropic.AsyncAnthropic(api_key=self.api_key)

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
