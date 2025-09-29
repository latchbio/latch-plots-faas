import asyncio
import os
import socket
import sys
import traceback
import uuid
from dataclasses import asdict, dataclass, field
from enum import Enum
from textwrap import dedent

from agents import Agent, Runner, SQLiteSession, function_tool
from agents.model_settings import ModelSettings
from instructions import construct_instructions
from lplots import _inject
from openai.types.shared.reasoning import Reasoning
from socketio_thread import SocketIoThread


class Mode(Enum):
    planning = "planning"
    executing = "executing"
    debugging = "debugging"


class PlanStatus(Enum):
    todo = "todo"
    in_progress = "in_progress"
    done = "done"


class PlanAction(Enum):
    add = "add"
    update = "update"
    complete = "complete"


@dataclass
class PlanItem:
    id: str
    description: str
    status: PlanStatus


@dataclass
class PlanDiff:
    action: PlanAction
    id: str
    description: str


@dataclass
class NotebookResponse:
    plan: list[PlanItem]
    plan_diff: list[PlanDiff]
    summary: list[str] | None = None
    questions: list[str] | None = None


@dataclass(kw_only=True)
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
            await self.send({"type": "agent_action", "action": action, "params": params, "tx_id": tx_id})
        except Exception as e:
            self.pending_operations.pop(tx_id, None)
            return {"status": "error", "error": f"Send failed: {e!s}"}

        try:
            return await asyncio.wait_for(response_future, timeout=5.0)
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
                return f"Created cell at position {position} (ID: {cell_id}, Title: {title})"
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

            params = {
                "position": position,
                "cell_type": "markdown",
                "source": code,
            }

            result = await self.atomic_operation("create_markdown_cell", params)
            if result.get("status") == "success":
                cell_id = result.get("cell_id", "unknown")
                return f"Created markdown cell at position {position} (ID: {cell_id})"
            return f"Failed to create cell: {result.get('error', 'Unknown error')}"

        @function_tool
        async def edit_cell(cell_id: str, new_code: str, auto_run: bool = True) -> str:
            """Replace the contents of an existing cell.

            Args:
                cell_id: The cell ID - can be either the user-visible number (0, 1, 2...) or internal ID
                new_code: The new source code for the cell
                auto_run: Whether to automatically run the cell after editing
            """
            params = {
                "cell_id": cell_id,
                "source": new_code,
                "auto_run": auto_run
            }

            result = await self.atomic_operation("edit_cell", params)
            if result.get("status") == "success":
                return f"Cell {cell_id} edited successfully"
            return f"Failed to edit cell: {result.get('error', 'Unknown error')}"

        @function_tool
        async def delete_cell(cell_id: str) -> str:
            """Remove a cell from the notebook."""
            params = {"cell_id": cell_id}

            result = await self.atomic_operation("delete_cell", params)
            if result.get("status") == "success":
                remaining = result.get("remaining_cells", [])
                cell_count = result.get("cell_count", 0)

                if remaining:
                    cell_list = ", ".join([f"{c['index']}: {c['cell_type']}" for c in remaining[:5]])
                    if len(remaining) > 5:
                        cell_list += f", ... ({len(remaining) - 5} more)"
                    return f"Cell {cell_id} deleted. {cell_count} cells remain: [{cell_list}]"
                return f"Cell {cell_id} deleted. No cells remain in notebook."
            return f"Failed to delete cell: {result.get('error', 'Unknown error')}"

        @function_tool
        async def run_cell(cell_id: str) -> str:
            """Execute a specific cell."""
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
            """Stop execution of a specific cell."""
            params = {"cell_id": cell_id}

            result = await self.atomic_operation("stop_cell", params)
            if result.get("status") == "success":
                self.executing_cells.discard(cell_id)
                return f"Stopped cell {cell_id}"
            return f"Failed to stop cell {cell_id}: {result.get('error', 'Unknown error')}"

        @function_tool
        async def delete_all_cells() -> str:
            """Delete all cells in the notebook efficiently."""
            params = {}

            result = await self.atomic_operation("clear_notebook", params)
            if result.get("status") == "success":
                return "Cleared all cells from the notebook"
            return f"Failed to clear notebook: {result.get('error', 'Unknown error')}"

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
                user_id = cell.get("user_id", cell.get("index", "?"))
                cell_type = cell.get("cell_type", "unknown")
                status = cell.get("status", "idle")
                source = cell.get("source", "")

                source_preview = source[:100] + "..." if len(source) > 100 else source
                source_preview = source_preview.replace("\n", " ")

                summary += f"\nCell {user_id} ({cell_type}, {status})"
                if source_preview:
                    summary += f": {source_preview}"

            return summary

        @function_tool
        async def send_plan_update(
            plan: list[dict],
            plan_diff: list[dict] | None = None,
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
            try:
                plan_items = [PlanItem(**item) for item in plan]
                plan_diff_items = [PlanDiff(**item) for item in (plan_diff or [])]
            except Exception as e:
                return f"Invalid plan data: {e}"

            plan_payload = {
                "plan": [asdict(item) for item in plan_items],
                "plan_diff": [asdict(item) for item in plan_diff_items],
            }

            result = await self.atomic_operation("plan_update", plan_payload)
            if result.get("status") == "success":
                return "Plan update delivered"
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
        print("[agent] Initializing")

        self.api_key = os.environ.get("OPENAI_API_KEY")

        if self.api_key:
            try:
                os.environ["OPENAI_API_KEY"] = self.api_key

                self.session = SQLiteSession()
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
                print("[agent] Initialization complete")
            except Exception as e:
                await self.send({
                    "type": "agent_error",
                    "error": f"Failed to initialize: {e!s}",
                    "fatal": True
                })
        else:
            await self.send({
                "type": "agent_error",
                "error": "OPENAI_API_KEY not set",
                "fatal": True
            })

    async def handle_query(self, msg: dict[str, object]) -> None:
        query = msg.get("query", "")

        print(f"[agent] Processing query: {query[:100]}...")

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

            response_msg = {
                "type": "agent_result",
                "status": "success",
                "response": response_content,
                "mode": self.mode.value
            }

            if structured_output is not None:
                response_msg["structured_output"] = asdict(structured_output)

            await self.send(response_msg)

        except Exception as e:
            print(f"[agent] Error processing query: {e}")
            await self.send({
                "type": "agent_result",
                "status": "error",
                "error": str(e)
            })

    def create_tracked_task(self, coro) -> asyncio.Task:
        task = asyncio.create_task(coro)
        self.active_tasks.add(task)
        task.add_done_callback(self.active_tasks.discard)
        return task

    async def fix_cell_error(self, cell_id: str, exception: str) -> None:
        print(f"[agent] Auto-fixing error in cell {cell_id}")

        try:
            result = await self.atomic_operation("get_context", {})
            source = ""
            if result.get("status") == "success":
                cells = result.get("context", {}).get("cells", [])
                for cell in cells:
                    if str(cell.get("id")) == str(cell_id) or str(cell.get("user_id")) == str(cell_id):
                        source = cell.get("source", "")
                        break

            if source == "":
                print(f"[agent] Could not find source for cell {cell_id}")
                return

            fix_query = dedent(f"""
                Cell {cell_id} failed with this error:

                ```python
                {source}
                ```

                Error:
                ```
                {exception}
                ```

                Analyze and fix this error by:
                1. Understanding what went wrong
                2. Either editing the cell to fix it OR creating a new cell before it with necessary imports/setup
                3. Re-running the cell after fixing

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
                if self.mode != previous_mode:
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
        print("[agent] Cancel request")

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

        print(f"[agent] Received message type: {msg_type}")

        if msg_type == "init":
            await self.handle_init(msg)
        elif msg_type == "agent_query":
            await self.handle_query(msg)
        elif msg_type == "agent_cancel":
            await self.handle_cancel(msg)
        elif msg_type == "agent_action_response":
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
        else:
            print(f"[agent] Unknown message type: {msg_type}")


async def main() -> None:
    global loop
    loop = asyncio.get_running_loop()

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
