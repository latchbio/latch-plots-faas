import asyncio
import json
import os
import socket
import sys
import traceback
import uuid
from dataclasses import dataclass, field
from enum import Enum
from textwrap import dedent
from typing import Callable, Literal

import anthropic
from anthropic.types import MessageParam, ToolParam, ToolUseBlock
from config_loader import build_full_instruction, get_custom_tools
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

# todo(tim): anthropic offfers @beta_tool decorator and tool_runner, but not sure
# these would make sense given our custom reqs like mode switching, cell exec tracking, etc.
class ToolDefinition:
    def __init__(self, name: str, description: str, input_schema: dict, handler: Callable):
        self.name = name
        self.description = description
        self.input_schema = input_schema
        self.handler = handler

    def to_anthropic_tool(self) -> ToolParam:
        return {
            "name": self.name,
            "description": self.description,
            "input_schema": self.input_schema,
        }


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


@dataclass
class AgentHarness:
    conn: SocketIoThread
    initialized: bool = False
    api_key: str | None = None
    agent: anthropic.AsyncAnthropic | None = None
    conversation_history: list[MessageParam] = field(default_factory=list)
    mode: Mode = Mode.planning
    pending_operations: dict[str, asyncio.Future] = field(default_factory=dict)
    executing_cells: set[str] = field(default_factory=set)
    tools: list[ToolDefinition] = field(default_factory=list)
    tool_map: dict[str, ToolDefinition] = field(default_factory=dict)
    active_tasks: set[asyncio.Task] = field(default_factory=set)
    error_fixes_in_progress: set[str] = field(default_factory=set)
    operation_counter: int = 0
    instructions_context: str = ""

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

    def register_tool(self, tool: ToolDefinition) -> None:
        """Register tool for API exposure and quick lookup."""
        self.tools.append(tool)
        self.tool_map[tool.name] = tool
    
    def _last_assistant_has_thinking(self, messages: list[MessageParam]) -> bool:
        for msg in reversed(messages):
            if msg.get("role") != "assistant":
                continue

            content = msg.get("content", [])
            if not (isinstance(content, list) and content):
                return False

            first_block = content[0]
            if isinstance(first_block, dict):
                block_type = first_block.get("type")
            else:
                block_type = getattr(first_block, "type", None)

            return block_type in ("thinking", "redacted_thinking")

        return False

    def init_tools(self) -> None:
        self.tools.clear()
        self.tool_map.clear()

        async def create_cell_handler(
            position: int,
            code: str,
            title: str,
            auto_run: bool = True,
        ) -> str:
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

        self.register_tool(ToolDefinition(
            name="create_cell",
            description="Create a new code cell at specified position.",
            input_schema={
                "type": "object",
                "properties": {
                    "position": {
                        "type": "integer",
                        "description": "The position to insert the cell at (0-indexed)"
                    },
                    "code": {
                        "type": "string",
                        "description": "The source code to put in the cell"
                    },
                    "title": {
                        "type": "string",
                        "description": "Descriptive title (<=6 words, Title Case) for the cell"
                    },
                    "auto_run": {
                        "type": "boolean",
                        "description": "Whether to automatically run the cell after creation",
                        "default": True
                    }
                },
                "required": ["position", "code", "title"]
            },
            handler=create_cell_handler
        ))

        async def create_markdown_cell_handler(position: int, code: str) -> str:
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

        self.register_tool(ToolDefinition(
            name="create_markdown_cell",
            description="Create a new markdown cell at specified position.",
            input_schema={
                "type": "object",
                "properties": {
                    "position": {
                        "type": "integer",
                        "description": "The position to insert the cell at (0-indexed)"
                    },
                    "code": {
                        "type": "string",
                        "description": "The source markdown to put in the cell"
                    }
                },
                "required": ["position", "code"]
            },
            handler=create_markdown_cell_handler
        ))

        async def edit_cell_handler(cell_id: str, new_code: str, auto_run: bool = True) -> str:
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

        self.register_tool(ToolDefinition(
            name="edit_cell",
            description="Replace the contents of an existing cell.",
            input_schema={
                "type": "object",
                "properties": {
                    "cell_id": {
                        "type": "string",
                        "description": "The exact cell_id string from get_notebook_context (e.g. 'cid:0@123:Map'), NOT the index number."
                    },
                    "new_code": {
                        "type": "string",
                        "description": "The new source code for the cell"
                    },
                    "auto_run": {
                        "type": "boolean",
                        "description": "Whether to automatically run the cell after editing",
                        "default": True
                    }
                },
                "required": ["cell_id", "new_code"]
            },
            handler=edit_cell_handler
        ))

        async def delete_cell_handler(cell_id: str) -> str:
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

        self.register_tool(ToolDefinition(
            name="delete_cell",
            description="Remove a cell from the notebook.",
            input_schema={
                "type": "object",
                "properties": {
                    "cell_id": {
                        "type": "string",
                        "description": "The exact cell_id string from get_notebook_context (e.g. 'cid:0@123:Map'), NOT the index number."
                    }
                },
                "required": ["cell_id"]
            },
            handler=delete_cell_handler
        ))

        async def run_cell_handler(cell_id: str) -> str:
            params = {"cell_id": cell_id}

            await self.send({
                "type": "agent_action",
                "action": "run_cell",
                "params": params,
            })
            self.executing_cells.add(cell_id)

            return f"Cell {cell_id} execution started"

        self.register_tool(ToolDefinition(
            name="run_cell",
            description="Execute a specific cell.",
            input_schema={
                "type": "object",
                "properties": {
                    "cell_id": {
                        "type": "string",
                        "description": "The exact cell_id string from get_notebook_context (e.g. 'cid:0@123:Map'), NOT the index number."
                    }
                },
                "required": ["cell_id"]
            },
            handler=run_cell_handler
        ))

        async def stop_cell_handler(cell_id: str) -> str:
            params = {"cell_id": cell_id}

            result = await self.atomic_operation("stop_cell", params)
            if result.get("status") == "success":
                self.executing_cells.discard(cell_id)
                return f"Stopped cell {cell_id}"
            return f"Failed to stop cell {cell_id}: {result.get('error', 'Unknown error')}"

        self.register_tool(ToolDefinition(
            name="stop_cell",
            description="Stop execution of a specific cell.",
            input_schema={
                "type": "object",
                "properties": {
                    "cell_id": {
                        "type": "string",
                        "description": "The exact cell_id string from get_notebook_context (e.g. 'cid:0@123:Map'), NOT the index number."
                    }
                },
                "required": ["cell_id"]
            },
            handler=stop_cell_handler
        ))

        async def delete_all_cells_handler() -> str:
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

        self.register_tool(ToolDefinition(
            name="delete_all_cells",
            description="Delete all cells in the notebook efficiently.",
            input_schema={
                "type": "object",
                "properties": {}
            },
            handler=delete_all_cells_handler
        ))

        async def get_notebook_context_handler() -> str:
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

                source_preview = source[:500] + "..." if len(source) > 500 else source
                source_preview = source_preview.replace("\n", " ")

                summary += f"\n[{index}] ({cell_type}, {status}, id: {cell_id})"
                if source_preview:
                    summary += f": {source_preview}"

            return summary

        self.register_tool(ToolDefinition(
            name="get_notebook_context",
            description="Get the current state of the notebook including all cells and their content.",
            input_schema={
                "type": "object",
                "properties": {}
            },
            handler=get_notebook_context_handler
        ))

        async def send_plan_update_handler(
            plan: list[PlanItemPayload],
            plan_diff: list[PlanDiffPayload],
        ) -> str:
            if AGENT_DEBUG:
                print(f"[tool] send_plan_update plan_items={len(plan)} diff_items={len(plan_diff)}")

            try:
                plan_items = [PlanItem(**item) for item in plan]
                plan_diff_items = [PlanDiff(**item) for item in plan_diff]
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

        self.register_tool(ToolDefinition(
            name="send_plan_update",
            description="Send the latest non-final plan state to the frontend. The final agent response already includes the complete plan, so skip this call if you are about to return the final response.",
            input_schema={
                "type": "object",
                "properties": {
                    "plan": {
                        "type": "array",
                        "description": "Full set of current plan items",
                        "items": {
                            "type": "object",
                            "properties": {
                                "id": {"type": "string"},
                                "description": {"type": "string"},
                                "status": {"type": "string", "enum": ["todo", "in_progress", "done"]}
                            },
                            "required": ["id", "description", "status"]
                        }
                    },
                    "plan_diff": {
                        "type": "array",
                        "description": "Updates describing what changed. Provide an empty array when nothing changed since the last update.",
                        "items": {
                            "type": "object",
                            "properties": {
                                "action": {"type": "string", "enum": ["add", "update", "complete"]},
                                "id": {"type": "string"},
                                "description": {"type": "string"}
                            },
                            "required": ["action", "id", "description"]
                        }
                    }
                },
                "required": ["plan", "plan_diff"]
            },
            handler=send_plan_update_handler
        ))

        async def start_new_plan_handler() -> str:
            self.set_mode(Mode.planning)
            return "Started new planning session"

        self.register_tool(ToolDefinition(
            name="start_new_plan",
            description="Start a new planning session. Only use this tool if the previous plan is complete and the user asks for tasks that require a new plan.",
            input_schema={
                "type": "object",
                "properties": {}
            },
            handler=start_new_plan_handler
        ))

        async def submit_response_handler(
            plan: list[PlanItemPayload],
            plan_diff: list[PlanDiffPayload],
            summary: list[str] | None = None,
            questions: list[str] | None = None,
        ) -> str:
            """This is a special tool that doesn't actually execute - it's used to capture structured output."""
            return "Response submitted"

        self.register_tool(ToolDefinition(
            name="submit_response",
            description="Submit your final response to the user with a structured plan, plan updates, optional summary points, and optional questions. Use this tool when you are ready to respond to the user.",
            input_schema={
                "type": "object",
                "properties": {
                    "plan": {
                        "type": "array",
                        "description": "Full set of current plan items showing all tasks",
                        "items": {
                            "type": "object",
                            "properties": {
                                "id": {"type": "string", "description": "Unique identifier for the plan item"},
                                "description": {"type": "string", "description": "Description of the task"},
                                "status": {"type": "string", "enum": ["todo", "in_progress", "done"], "description": "Current status of the task"}
                            },
                            "required": ["id", "description", "status"]
                        }
                    },
                    "plan_diff": {
                        "type": "array",
                        "description": "Changes to the plan since last update. Use an empty array when nothing changed.",
                        "items": {
                            "type": "object",
                            "properties": {
                                "action": {"type": "string", "enum": ["add", "update", "complete"], "description": "Type of change"},
                                "id": {"type": "string", "description": "ID of the affected plan item"},
                                "description": {"type": "string", "description": "Description of the change"}
                            },
                            "required": ["action", "id", "description"]
                        }
                    },
                    "summary": {
                        "type": "array",
                        "description": "Optional summary points about what was done or discussed",
                        "items": {"type": "string"}
                    },
                    "questions": {
                        "type": "array",
                        "description": "Optional questions for the user",
                        "items": {"type": "string"}
                    }
                },
                "required": ["plan", "plan_diff"]
            },
            handler=submit_response_handler
        ))

    def _get_system_prompt(self) -> str:
        return build_full_instruction(self.instructions_context)

    async def _call_tool(self, tool_name: str, tool_input: dict) -> str:
        tool = self.tool_map.get(tool_name)
        if tool is None:
            return f"Unknown tool: {tool_name}"

        try:
            result = await tool.handler(**tool_input)
            return str(result)
        except Exception as e:
            error_msg = f"Error executing {tool_name}: {str(e)}"
            print(f"[agent] {error_msg}")
            if AGENT_DEBUG:
                traceback.print_exc()
            return error_msg

    async def _run_agent_step(self, user_query: str) -> tuple[list[str], NotebookResponse | None]:
        responses = []

        self.conversation_history.append({
            "role": "user",
            "content": user_query
        })

        model, thinking_budget = self.mode_config[self.mode]
        max_turns = 100
        turn = 0

        while turn < max_turns:
            turn += 1

            try:
                if thinking_budget is not None:
                    max_tokens = thinking_budget + 4096
                else:
                    max_tokens = 4096

                # Filter conversation history when thinking is enabled
                messages_to_send = self.conversation_history
                if thinking_budget is not None:
                    valid_tool_use_ids = set()
                    filtered_messages = []
                    removed_count = 0
                    filtered_assistant_count = 0
                    assistant_messages_seen = any(
                        msg.get("role") == "assistant" for msg in self.conversation_history
                    )

                    def _block_type(block):
                        if isinstance(block, dict):
                            return block.get("type")
                        return getattr(block, "type", None)

                    def _tool_use_id(block):
                        if isinstance(block, dict):
                            return block.get("id")
                        return getattr(block, "id", None)

                    # todo(tim): clean this up, try a one pass approach instead
                    for msg in self.conversation_history:
                        if msg["role"] == "assistant":
                            content = msg.get("content", [])
                            if isinstance(content, list) and len(content) > 0:
                                has_tool_use = False
                                for block in content:
                                    if _block_type(block) == "tool_use":
                                        has_tool_use = True
                                        tool_id = _tool_use_id(block)
                                        if tool_id:
                                            valid_tool_use_ids.add(tool_id)

                    for msg in self.conversation_history:
                        if msg["role"] == "assistant":
                            content = msg.get("content", [])
                            if isinstance(content, list) and len(content) > 0:
                                first_block = content[0]
                                has_thinking = (
                                    _block_type(first_block) in ("thinking", "redacted_thinking")
                                )

                                has_tool_use = False
                                for block in content:
                                    if _block_type(block) == "tool_use":
                                        has_tool_use = True
                                        break

                                if has_thinking or has_tool_use:
                                    filtered_messages.append(msg)
                                    filtered_assistant_count += 1
                                else:
                                    removed_count += 1
                        elif msg["role"] == "user":
                            content = msg.get("content", [])
                            if isinstance(content, list):
                                filtered_content = []
                                for block in content:
                                    if isinstance(block, dict) and block.get("type") == "tool_result":
                                        tool_use_id = block.get("tool_use_id")
                                        if tool_use_id in valid_tool_use_ids:
                                            filtered_content.append(block)
                                    else:
                                        filtered_content.append(block)

                                if filtered_content:
                                    filtered_messages.append({
                                        "role": "user",
                                        "content": filtered_content
                                    })
                            else:
                                filtered_messages.append(msg)

                    messages_to_send = filtered_messages
                    if filtered_assistant_count == 0 and assistant_messages_seen:
                        messages_to_send = self.conversation_history
                        removed_count = 0

                request_params = {
                    "model": model,
                    "max_tokens": max_tokens,
                    "system": self._get_system_prompt(),
                    "messages": messages_to_send,
                    "tools": [tool.to_anthropic_tool() for tool in self.tools],
                }

                enable_thinking = thinking_budget is not None
                if enable_thinking:
                    request_params["thinking"] = {
                        "type": "enabled",
                        "budget_tokens": thinking_budget
                    }
                    # note(tim): anthropic requires thinking blocks at the beginning of the assistant message
                    if not self._last_assistant_has_thinking(messages_to_send):
                        enable_thinking = False
                        request_params.pop("thinking", None)

                response = await self.agent.messages.create(**request_params)


            except Exception as e:
                raise RuntimeError(f"Error calling Claude API: {e}")

            # Add assistant response to history
            assistant_content = []
            tool_uses = []
            text_responses = []

            for block in response.content:
                if block.type == "text":
                    text_responses.append(block.text)
                    assistant_content.append(block)
                elif block.type == "tool_use":
                    tool_uses.append(block)
                    assistant_content.append(block)

            self.conversation_history.append({
                "role": "assistant",
                "content": assistant_content
            })

            if not tool_uses:
                if text_responses:
                    combined_text = "\n".join(text_responses)
                    responses.append(combined_text)

                return responses, None

            tool_results = []
            structured_output_from_tool = None

            for tool_use in tool_uses:
                tool_name = tool_use.name
                tool_input = tool_use.input

                if tool_name == "submit_response":
                    print(f"[agent] Received submit_response with structured output")
                    try:
                        summary = tool_input.get("summary")
                        if summary is not None:
                            if isinstance(summary, str) and summary.lower() in ("null", "none"):
                                summary = None
                            elif not isinstance(summary, list):
                                print(f"[agent] Warning: summary is not a list, got {type(summary)}")
                                summary = None

                        questions = tool_input.get("questions")
                        if questions is not None:
                            if isinstance(questions, str) and questions.lower() in ("null", "none"):
                                questions = None
                            elif not isinstance(questions, list):
                                print(f"[agent] Warning: questions is not a list, got {type(questions)}")
                                questions = None

                        structured_output_from_tool = NotebookResponse(
                            plan=[PlanItem(**item) for item in tool_input.get("plan", [])],
                            plan_diff=[PlanDiff(**item) for item in tool_input.get("plan_diff", [])],
                            summary=summary,
                            questions=questions
                        )
                        result = await self._call_tool(tool_name, tool_input)
                    except Exception as e:
                        print(f"[agent] Error parsing submit_response: {e}")
                        if AGENT_DEBUG:
                            traceback.print_exc()
                        result = f"Error parsing response: {e}"
                else:
                    if AGENT_DEBUG:
                        print(f"[agent] Executing tool: {tool_name}")
                    result = await self._call_tool(tool_name, tool_input)

                tool_results.append({
                    "type": "tool_result",
                    "tool_use_id": tool_use.id,
                    "content": result
                })

            self.conversation_history.append({
                "role": "user",
                "content": tool_results
            })

            # If we got structured output from submit_response, we're done
            if structured_output_from_tool is not None:
                if text_responses:
                    responses.append("\n".join(text_responses))

                return responses, structured_output_from_tool

        # Max turns reached
        responses.append("Maximum conversation turns reached without completion")
        return responses, None

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
            self.agent = anthropic.AsyncAnthropic(api_key=self.api_key)
            self.init_tools()

            context = msg.get("context", "")
            self.instructions_context = context

            # Get initial notebook context
            try:
                result = await self.atomic_operation("get_context", {})
                if result.get("status") == "success":
                    ctx = result.get("context", {})
                    cell_count = ctx.get("cell_count", 0)
                    cells = ctx.get("cells", [])

                    if cell_count > 0:
                        self.instructions_context += f"\n\nCurrent notebook state: {cell_count} cell(s) present."
                        for cell in cells[:5]:
                            idx = cell.get("index", "?")
                            cell_type = cell.get("cell_type", "unknown")
                            source = cell.get("source", "")[:50]
                            if source:
                                self.instructions_context += f"\n- Cell {idx} ({cell_type}): {source}..."
                        if cell_count > 5:
                            self.instructions_context += f"\n- ... and {cell_count - 5} more cells"
            except Exception as e:
                print(f"[agent] Warning: Could not fetch initial context: {e}")

            self.initialized = True
            await self.send({
                "type": "agent_status",
                "status": "ready"
            })
            print("[agent] Initialization complete", flush=True)
        except Exception as e:
            await self.send({
                "type": "agent_error",
                "error": f"Failed to initialize: {e!s}",
                "fatal": True
            })

    async def handle_query(self, msg: dict[str, object]) -> None:
        query = msg.get("query", "")
        request_id = msg.get("request_id")

        print(f"[agent] Processing query: {query[:500]}...")


        try:
            previous_ops = self.operation_counter

            if self.mode == Mode.planning:
                responses, structured_output = await self._run_agent_step(query)
                if self.operation_counter > previous_ops:
                    self.set_mode(Mode.executing)
            else:
                if self.mode != Mode.executing:
                    self.set_mode(Mode.executing)
                responses, structured_output = await self._run_agent_step(query)

            if self.mode == Mode.executing:
                self.set_mode(Mode.planning)

            response_msg = {
                "type": "agent_result",
                "status": "success",
                "responses": responses if responses else [""],
                "mode": self.mode.value,
            }

            if request_id is not None:
                response_msg["request_id"] = request_id


            structured_output_dict = None

            if structured_output is not None:
                structured_output_dict = structured_output.model_dump()
                # TODO(tim): cleanup this code, consider deleting
                # Ensure required plan fields are always present as arrays
                if not isinstance(structured_output_dict.get("plan"), list):
                    structured_output_dict["plan"] = []
                if not isinstance(structured_output_dict.get("plan_diff"), list):
                    structured_output_dict["plan_diff"] = []
            # todo(tim): cleanup this code, consider deleting
            else:
                structured_output_dict = {
                    "plan": [],
                    "plan_diff": [],
                    "summary": None,
                    "questions": None,
                }
            # todo(tim): cleanup this code, consider deleting just his conditional (not indisde it)
            if structured_output_dict is not None:
                # Ensure summary and questions are None or arrays
                if structured_output_dict.get("summary") is not None and not isinstance(structured_output_dict["summary"], list):
                    structured_output_dict["summary"] = None
                if structured_output_dict.get("questions") is not None and not isinstance(structured_output_dict["questions"], list):
                    structured_output_dict["questions"] = None

                # If there's no summary but there are text responses, put them in summary
                if (structured_output_dict.get("summary") is None or len(structured_output_dict.get("summary", [])) == 0):
                    if responses and len(responses) > 0:
                        structured_output_dict["summary"] = responses

                response_msg["structured_output"] = structured_output_dict
            print(f"[agent] Response message: {response_msg}")

            await self.send(response_msg)

        except Exception as e:
            print(f"[agent] Error processing query: {e}")
            if AGENT_DEBUG:
                traceback.print_exc()
            error_msg = {
                "type": "agent_result",
                "status": "error",
                "error": str(e),
                "mode": self.mode.value,
                "structured_output": {
                    "plan": [],
                    "plan_diff": [],
                    "summary": None,
                    "questions": None,
                },
            }

            if request_id is not None:
                error_msg["request_id"] = request_id

            await self.send(error_msg)

    def create_tracked_task(self, coro) -> asyncio.Task:
        task = asyncio.create_task(coro)
        self.active_tasks.add(task)
        task.add_done_callback(self.active_tasks.discard)
        return task

    async def fix_cell_error(self, cell_id: str, exception: str) -> None:
        print(f"[agent] Auto-fixing error in cell {cell_id}")

        if not exception:
            print(f"[agent] No exception text provided for cell {cell_id}")
            return

        try:
            ctx_result = await self.atomic_operation("get_context", {})

            if ctx_result.get("status") != "success":
                print("[agent] Could not fetch context")
                return

            cells = ctx_result.get("context", {}).get("cells", [])
            cell_info = next((c for c in cells if str(c.get("tf_id")) == str(cell_id)), None)

            if cell_info is None:
                print(f"[agent] Cell with tf_id={cell_id} not found in context")
                print(f"[agent] Available cells: {[(c.get('cell_id'), c.get('tf_id')) for c in cells]}")
                return

            source = cell_info.get("source", "")
            crdt_cell_id = cell_info.get("cell_id")

            if source == "":
                print(f"[agent] Cell {cell_id} has no source code")
                return

            if crdt_cell_id is None:
                print(f"[agent] No cell_id for tf_id {cell_id}")
                return

            exception_text = exception
            try:
                parsed = json.loads(exception)
                if isinstance(parsed, dict) and "string" in parsed:
                    exception_text = parsed["string"]
            except json.JSONDecodeError:
                print("[agent] Exception is not JSON-encoded, using raw text")

            fix_query = dedent(f"""
                Cell {crdt_cell_id} at position {cell_info.get('index', '?')} failed with this error:

                ```python
                {source}
                ```

                Error:
                ```
                {exception_text}
                ```

                Fix this error by editing the cell with the corrected code.
                Only create additional cells if the error is due to missing
                imports or dependencies that should be in a separate cell.

                Be concise and fix the issue directly.
            """)

            previous_mode = self.mode
            self.set_mode(Mode.debugging)

            try:
                responses, _ = await self._run_agent_step(fix_query)

                await self.send({
                    "type": "agent_error_fixed",
                    "cell_id": cell_id,
                    "status": "success"
                })
                print(f"[agent] Successfully fixed error in cell {cell_id}")

            finally:
                self.set_mode(previous_mode)

        except Exception as e:
            print(f"[agent] Failed to fix error in cell {cell_id}: {e}")
            await self.send({
                "type": "agent_error_fixed",
                "cell_id": cell_id,
                "status": "failed",
                "error": str(e)
            })
        finally:
            self.error_fixes_in_progress.discard(cell_id)

    async def handle_cell_error(self, msg: dict[str, object]) -> None:
        cell_id = str(msg.get("cell_id"))
        exception = msg.get("exception", "")

        if cell_id in self.error_fixes_in_progress:
            print(f"[agent] Error fix for cell {cell_id} already in progress")
            return

        self.error_fixes_in_progress.add(cell_id)
        self.create_tracked_task(self.fix_cell_error(cell_id, exception))
        print(f"[agent] Started error fix task for cell {cell_id}")

    async def handle_cancel(self, msg: dict[str, object]) -> None:
        request_id = msg.get("request_id", "unknown")
        print(f"[agent] Cancelling request {request_id}")

        for task in list(self.active_tasks):
            if not task.done():
                task.cancel()

        if self.active_tasks:
            await asyncio.gather(*self.active_tasks, return_exceptions=True)

        self.active_tasks.clear()
        self.error_fixes_in_progress.clear()

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
            self.create_tracked_task(self.handle_query(msg))
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
        elif msg_type == "cell_result" and msg.get("exception") is not None:
            await self.handle_cell_error({
                "cell_id": msg.get("cell_id"),
                "exception": msg.get("exception")
            })
        elif msg_type == "start_cell":
            cell_id = msg.get("cell_id")
            if cell_id is not None:
                self.executing_cells.add(str(cell_id))
        elif msg_type == "cell_result":
            cell_id = msg.get("cell_id")
            if cell_id is not None:
                self.executing_cells.discard(str(cell_id))
        elif msg_type == "kernel_message":
            nested_msg = msg.get("message", {})
            nested_type = nested_msg.get("type")
            if nested_type == "cell_result" and nested_msg.get("has_exception") is True:
                cell_id = nested_msg.get("cell_id")
                exception = nested_msg.get("exception", "")
                print(f"[agent] Cell error detected: cell_id={cell_id}, has_exception={exception != ''}")
                await self.handle_cell_error({
                    "cell_id": cell_id,
                    "exception": exception
                })
            elif nested_type == "start_cell":
                cell_id = nested_msg.get("cell_id")
                if cell_id is not None:
                    self.executing_cells.add(str(cell_id))
            elif nested_type == "cell_result":
                cell_id = nested_msg.get("cell_id")
                if cell_id is not None:
                    self.executing_cells.discard(str(cell_id))
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
