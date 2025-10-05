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
from typing import Literal

from agents import Agent, Runner, SQLiteSession, function_tool
from agents.model_settings import ModelSettings
from instructions import construct_instructions
from lplots import _inject
from openai.types.shared.reasoning import Reasoning
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


@dataclass
class AgentHarness:
    conn: SocketIoThread
    initialized: bool = False
    api_key: str | None = None
    agent: Agent | None = None
    session: SQLiteSession | None = None
    mode: Mode = Mode.planning
    pending_operations: dict[str, asyncio.Future] = field(default_factory=dict)
    executing_cells: set[str] = field(default_factory=set)
    tools: list = field(default_factory=list)
    active_tasks: set[asyncio.Task] = field(default_factory=set)
    error_fixes_in_progress: set[str] = field(default_factory=set)
    operation_counter: int = 0

    mode_config: dict[Mode, tuple[str, str]] = field(default_factory=lambda: {
        Mode.planning: ("gpt-5", "medium"),
        Mode.executing: ("gpt-5", "low"),
        Mode.debugging: ("gpt-5", "medium"),
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

        assert self.agent is not None, "Agent not initialized"

        model, reasoning_effort = self.mode_config[mode]
        self.mode = mode

        self.agent.model = model
        self.agent.model_settings = ModelSettings(
            reasoning=Reasoning(effort=reasoning_effort)
        )
        print(f"[agent] Mode changed to {mode.value}")

    def init_tools(self) -> None:
        @function_tool
        async def create_cell(
            position: int,
            code: str,
            title: str,
            auto_run: bool = True,
        ) -> str:
            """Create a new code cell at specified position.

            Args:
                position: The position to insert the cell at (0-indexed)
                code: The source code to put in the cell
                title: Descriptive title (<=6 words, Title Case) for the cell
                auto_run: Whether to automatically run the cell after creation
            """
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

        @function_tool
        async def create_markdown_cell(position: int, code: str) -> str:
            """Create a new markdown cell at specified position.

            Args:
                position: The position to insert the cell at (0-indexed)
                code: The source markdown to put in the cell
            """
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

        @function_tool
        async def edit_cell(cell_id: str, new_code: str, auto_run: bool = True) -> str:
            """Replace the contents of an existing cell.

            Args:
                cell_id: The exact cell_id string from get_notebook_context (e.g. 'cid:0@123:Map'), NOT the index number.
                new_code: The new source code for the cell
                auto_run: Whether to automatically run the cell after editing
            """
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

        @function_tool
        async def delete_cell(cell_id: str) -> str:
            """Remove a cell from the notebook.

            Args:
                cell_id: The exact cell_id string from get_notebook_context (e.g. 'cid:0@123:Map'), NOT the index number.
            """
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

        @function_tool
        async def run_cell(cell_id: str) -> str:
            """Execute a specific cell.

            Args:
                cell_id: The exact cell_id string from get_notebook_context (e.g. 'cid:0@123:Map'), NOT the index number.
            """
            params = {"cell_id": cell_id}

            await self.send({
                "type": "agent_action",
                "action": "run_cell",
                "params": params,
            })
            self.executing_cells.add(cell_id)

            return f"Cell {cell_id} execution started"

        @function_tool
        async def stop_cell(cell_id: str) -> str:
            """Stop execution of a specific cell.

            Args:
                cell_id: The exact cell_id string from get_notebook_context (e.g. 'cid:0@123:Map'), NOT the index number.
            """
            params = {"cell_id": cell_id}

            result = await self.atomic_operation("stop_cell", params)
            if result.get("status") == "success":
                self.executing_cells.discard(cell_id)
                return f"Stopped cell {cell_id}"
            return f"Failed to stop cell {cell_id}: {result.get('error', 'Unknown error')}"

        @function_tool
        async def delete_all_cells() -> str:
            """Delete all cells in the notebook efficiently."""
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

        @function_tool
        async def get_notebook_context() -> str:
            """Get the current state of the notebook including all cells and their content."""
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

        @function_tool
        async def send_plan_update(
            plan: list[PlanItemPayload],
            plan_diff: list[PlanDiffPayload] | None = None,
        ) -> str:
            """Send plan state update to the frontend.

            Args:
                plan: Full set of current plan items
                plan_diff: Optional updates describing what changed

            Returns:
                A plan status message

            Note:
                The final agent response already includes the complete plan, so
                skip this call if you are about to return the final response.
            """
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

        @function_tool
        async def start_new_plan() -> str:
            """Start a new planning session."""
            self.set_mode(Mode.planning)
            return "Started new planning session"

        self.tools = [
            create_cell,
            create_markdown_cell,
            edit_cell,
            delete_cell,
            run_cell,
            stop_cell,
            delete_all_cells,
            get_notebook_context,
            send_plan_update,
            start_new_plan,
        ]

    async def handle_init(self, msg: dict[str, object]) -> None:
        print("[agent] Initializing", flush=True)

        self.api_key = os.environ.get("OPENAI_API_KEY")

        if not self.api_key:
            await self.send({
                "type": "agent_error",
                "error": "OPENAI_API_KEY not set",
                "fatal": True
            })
            return

        try:
            os.environ["OPENAI_API_KEY"] = self.api_key

            session_id = f"local_session_{uuid.uuid4().hex[:8]}"
            self.session = SQLiteSession(session_id)
            self.init_tools()

            context = msg.get("context", "")
            model, reasoning_effort = self.mode_config[self.mode]

            self.agent = Agent(
                name="NotebookAssistant",
                model=model,
                instructions=construct_instructions(context),
                tools=self.tools,
                output_type=NotebookResponse,
                model_settings=ModelSettings(
                    reasoning=Reasoning(effort=reasoning_effort)
                )
            )

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

        assert self.agent is not None, "Agent not initialized"

        try:
            previous_ops = self.operation_counter

            result = await Runner.run(
                self.agent,
                query,
                session=self.session,
                max_turns=100
            )

            if self.mode == Mode.planning and self.operation_counter > previous_ops:
                self.set_mode(Mode.executing)

            if AGENT_DEBUG and hasattr(result, "new_items"):
                for item in result.new_items:
                    if hasattr(item, "raw_item") and hasattr(item.raw_item, "reasoning"):
                        reasoning = item.raw_item.reasoning
                        if reasoning and hasattr(reasoning, "content"):
                            print(f"[reasoning] {reasoning.content}")

            response_content = ""
            structured_output = None

            if hasattr(result, "final_output_as"):
                try:
                    structured_output = result.final_output_as(NotebookResponse)
                except Exception as e:
                    print(f"[agent] Could not extract structured output: {e}")

            if hasattr(result, "content"):
                response_content = str(result.content)
            elif hasattr(result, "output"):
                response_content = str(result.output)
            else:
                response_content = str(result)

            if self.mode == Mode.executing:
                self.set_mode(Mode.planning)

            response_msg = {
                "type": "agent_result",
                "status": "success",
                "responses": [response_content],
                "mode": self.mode.value,
            }

            if request_id is not None:
                response_msg["request_id"] = request_id

            if structured_output is not None:
                response_msg["structured_output"] = structured_output.model_dump()

            await self.send(response_msg)

        except Exception as e:
            print(f"[agent] Error processing query: {e}")
            await self.send({
                "type": "agent_result",
                "status": "error",
                "error": str(e),
                "mode": self.mode.value,
            })

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
                result = await Runner.run(
                    self.agent,
                    fix_query,
                    session=self.session,
                    max_turns=10
                )

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
