import asyncio
import json
import os
import socket
import subprocess  # noqa: S404
import sys
import time
import traceback
import uuid
from collections.abc import Callable, Iterable
from contextlib import suppress
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any

try:
    from claude_agent_sdk import ClaudeAgentOptions, ClaudeSDKClient
    from claude_agent_sdk.types import (
        AssistantMessage,
        ResultMessage,
        StreamEvent,
        SystemMessage,
        TextBlock,
        ThinkingBlock,
        ToolUseBlock,
    )
except ImportError:
    print("[agent] claude_agent_sdk missing, installing...", flush=True)
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", "claude-agent-sdk"])
    from claude_agent_sdk import ClaudeAgentOptions, ClaudeSDKClient
    from claude_agent_sdk.types import (
        AssistantMessage,
        ResultMessage,
        StreamEvent,
        SystemMessage,
        TextBlock,
        ThinkingBlock,
        ToolUseBlock,
    )
from lplots import _inject
from socketio_thread import SocketIoThread
from utils import auth_token_sdk, nucleus_url, sdk_token

sys.stdout.reconfigure(line_buffering=True)

skip_db_history = os.environ.get("AGENT_SKIP_DB_HISTORY") == "1"

cache_chunk_size = 20

reactivity_ready_statuses = {"ran", "ok", "success", "error"}

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


class Behavior(Enum):
    proactive = "proactive"
    step_by_step = "step_by_step"


@dataclass
class AgentHarness:
    conn: SocketIoThread
    initialized: bool = False
    client: ClaudeSDKClient | None = None
    mode: Mode = Mode.planning
    system_prompt: str | None = None
    pending_operations: dict[str, asyncio.Future] = field(default_factory=dict)
    executing_cells: set[str] = field(default_factory=set)
    tools: list[dict[str, object]] = field(default_factory=list)
    tool_map: dict[str, Callable] = field(default_factory=dict)
    operation_counter: int = 0
    current_request_id: str | None = None
    should_auto_continue: bool = False
    pending_auto_continue: bool = False
    pause_until_user_query: bool = False
    pending_tool_calls: set[str] = field(default_factory=set)

    pending_messages: asyncio.Queue = field(default_factory=asyncio.Queue)
    conversation_task: asyncio.Task | None = None
    conversation_running: bool = False
    agent_session_id: str | None = None
    latest_notebook_context: dict = field(default_factory=dict)
    current_status: str | None = None
    expected_widgets: dict[str, object | None] = field(default_factory=dict)
    behavior: Behavior = Behavior.step_by_step
    latest_notebook_state: str | None = None
    current_plan: dict | None = None
    buffer: list[str] = field(default_factory=list)
    summarize_tasks: set[asyncio.Task] = field(default_factory=set)
    in_memory_history: list[dict] = field(default_factory=list)
    mcp_servers: dict[str, object] = field(default_factory=dict)
    current_query_task: asyncio.Task | None = None
    open_stream_blocks: dict[int, str] = field(default_factory=dict)

    mode_config: dict[Mode, tuple[str, int | None]] = field(
        default_factory=lambda: {
            Mode.planning: ("claude-opus-4-5-20251101", 4096),
            Mode.executing: ("claude-opus-4-5-20251101", 1024),
            Mode.debugging: ("claude-opus-4-5-20251101", 2048),
        }
    )

    async def send(self, msg: dict[str, object]) -> None:
        msg_type = msg.get("type", "unknown")
        if msg_type == "agent_error":
            print(f"[agent] Sending message: agent_error payload={json.dumps(msg)}")
        elif msg_type != "agent_stream_delta":
            print(f"[agent] Sending message: {msg_type}")

        await self.conn.send(msg)

    async def refresh_cells_context(self) -> str:
        context_result, reactivity_result = await asyncio.gather(
            self.atomic_operation("get_context"),
            self.atomic_operation("request_reactivity_summary", timeout=1.0),
        )

        if context_result.get("status") != "success":
            print(f"[agent] Failed to get notebook context: {context_result}")
            return "Latch Plots is unable to provide context for this notebook due to an unknown error. Please inform the user that Latch Plots is having an issue and they should report it to support."

        reactivity_available: bool = reactivity_result.get("status") == "success"
        if not reactivity_available:
            print(
                f"[agent] Reactivity summary unavailable: {reactivity_result.get('error', 'unknown')}"
            )

        context = context_result.get("context", {})
        self.latest_notebook_context = context

        notebook_name = context.get("notebook_name")
        cell_count = context.get("cell_count", 0)
        cells = context.get("cells", [])

        cell_lines = [
            f"# Notebook Cells for {notebook_name}, Total cells: {cell_count}\n"
        ]

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

            if cell_type == "code":
                cell_name = cell.get("display_name", "Unknown Code Cell")
            else:
                cell_name = cell_type

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

            cell_lines.append(f"\n## Cell {cell_name} [{index}]")
            if cell_type == "code":
                cell_lines.append(f"DISPLAY_NAME: {cell_name}")
            cell_lines.append(f"BELONGS_TO_TAB: {current_tab_name}")  # noqa: FURB113
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

            if tf_id is None:
                continue

            cell_lines.append("\nREACTIVITY:")

            if not reactivity_available:
                cell_lines.append("- Reactivity summary unavailable.")
                continue

            cell_reactivity: dict[str, dict[str, Iterable[str]]] = (
                reactivity_result.get("cell_reactivity", {})
            )
            reactivity_meta = cell_reactivity.get(str(tf_id))
            is_reactivity_ready = status in reactivity_ready_statuses

            if not is_reactivity_ready:
                cell_lines.append("- Reactivity summary not ready.")
                continue

            if reactivity_meta is None:
                cell_lines.append("- Reactivity data missing for this cell.")
                continue

            signals_defined = reactivity_meta.get("signals_defined", [])
            depends_on_signals = reactivity_meta.get("depends_on_signals", [])
            depends_on_cells = reactivity_meta.get("depends_on_cells", [])

            cell_lines.append(  # noqa: FURB113
                "- Signals defined: "
                + (", ".join(signals_defined) if signals_defined else "None")
            )
            cell_lines.append(
                "- Depends on signals: "
                + (", ".join(depends_on_signals) if depends_on_signals else "None")
            )
            cell_lines.append(
                "- Depends on cells: "
                + (", ".join(depends_on_cells) if depends_on_cells else "None")
            )

        return "\n".join(cell_lines)

    def _behavior_file_name(self) -> str:
        return (
            "proactive.md"
            if self.behavior == Behavior.proactive
            else "step_by_step.md"
        )

    def _compose_turn_system_prompt(self) -> str:
        self.system_prompt = (context_root.parent / "system_prompt.md").read_text()
        assert self.system_prompt is not None

        behavior_file = self._behavior_file_name()
        turn_behavior_content = (context_root / "turn_behavior" / behavior_file).read_text()
        examples_content = (context_root / "examples" / behavior_file).read_text()

        final_system_prompt = self.system_prompt.replace(
            "TURN_BEHAVIOR_PLACEHOLDER",
            f"<turn_behavior>\n{turn_behavior_content}\n</turn_behavior>",
        )
        return final_system_prompt.replace(
            "EXAMPLES_PLACEHOLDER",
            f"<examples>\n{examples_content}\n</examples>",
        )

    async def _build_turn_prompt(self, user_query: str) -> str:
        notebook_state = await self.refresh_cells_context()
        self.latest_notebook_state = notebook_state

        context_blocks = [
            f"<current_notebook_state>\n{notebook_state}\n</current_notebook_state>"
        ]
        if self.current_plan is not None:
            plan_content = json.dumps(self.current_plan, indent=2)
            context_blocks.append(f"<current_plan>\n{plan_content}\n</current_plan>")

        context_blocks.append(f"<user_request>\n{user_query}\n</user_request>")
        return "\n\n".join(context_blocks)

    async def atomic_operation(
        self, action: str, params: dict | None = None, timeout: float = 10.0
    ) -> dict:
        if params is None:
            params: dict = {}

        self.operation_counter += 1

        if self.mode == Mode.planning and action in {
            "create_cell",
            "edit_cell",
            "run_cell",
            "delete_cell",
        }:
            self.set_mode(Mode.executing)

        force_backend_browser_retry = params.get("force_backend_browser_retry", False)

        tx_id = f"tx_{uuid.uuid4().hex[:12]}"
        loop = asyncio.get_running_loop()
        response_future = loop.create_future()
        self.pending_operations[tx_id] = response_future

        start_time = time.time()
        try:
            print(f"[agent] -> {action}")
            await self.send({
                "type": "agent_action",
                "action": action,
                "params": params,
                "tx_id": tx_id,
            })
        except Exception as e:
            self.pending_operations.pop(tx_id, None)
            return {"status": "error", "error": f"Send failed: {e!s}"}

        try:
            return await asyncio.wait_for(response_future, timeout=timeout)
        except asyncio.CancelledError:
            print(
                f"[agent] Operation cancelled (session reinitialized): action={action}, tx_id={tx_id}"
            )
            return {
                "status": "error",
                "error": f"OPERATION FAILED: '{action}' was interrupted because the session was reinitialized. This operation did NOT complete. You must retry or inform the user.",
            }
        except TimeoutError:
            self.pending_operations.pop(tx_id, None)
            duration = time.time() - start_time
            print(f"[agent] {action} timed out after {duration:.3f}s")

            if not force_backend_browser_retry:
                print(f"[agent] Retrying {action} with backend browser")
                return await self.atomic_operation(
                    action,
                    {**params, "force_backend_browser_retry": True},
                    timeout=timeout,
                )

            return {
                "status": "error",
                "error": f"OPERATION FAILED: '{action}' timed out after 10 seconds. This operation did NOT complete.",
                "tx_id": tx_id,
            }
        finally:
            ret = self.pending_operations.pop(tx_id, None)

            if ret is not None:
                duration = time.time() - start_time
                print(f"[agent] {action} took {duration:.3f}s")

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

    def init_tools(self) -> None:
        from tools import notebook_mcp

        self.mcp_servers = {"notebook-tools": notebook_mcp}

    async def _send_usage_update(self, usage: dict[str, Any]) -> None:
        await self.send({
            "type": "agent_usage_update",
            "input_tokens": int(usage.get("input_tokens", 0) or 0),
            "cache_read_input_tokens": int(usage.get("cache_read_input_tokens", 0) or 0),
            "cache_creation_input_tokens": int(
                usage.get("cache_creation_input_tokens", 0) or 0
            ),
            "context_limit": 200_000,
        })

    def _current_sdk_session_id(self) -> str:
        return self.agent_session_id or "default"

    @staticmethod
    def _normalize_session_id(raw_session_id: object) -> str | None:
        if raw_session_id is None:
            return None
        session_id = str(raw_session_id).strip()
        if session_id == "":
            return None
        return session_id

    async def _close_open_stream_blocks(self) -> None:
        if len(self.open_stream_blocks) == 0:
            return

        pending_blocks = sorted(self.open_stream_blocks)
        print(f"[agent] Closing {len(pending_blocks)} dangling stream block(s): {pending_blocks}")
        for block_index in pending_blocks:
            await self.send({"type": "agent_stream_block_stop", "block_index": block_index})
        self.open_stream_blocks.clear()

    async def _wait_for_running_query_to_stop(self, timeout_seconds: float = 2.0) -> bool:
        if self.current_query_task is None or self.current_query_task.done():
            return True

        try:
            await asyncio.wait_for(self.current_query_task, timeout=timeout_seconds)
            return True
        except TimeoutError:
            self.current_query_task.cancel()
            with suppress(asyncio.CancelledError):
                await self.current_query_task
            return False
        except asyncio.CancelledError:
            return True

    async def _reset_for_new_session(self) -> None:
        if len(self.pending_operations) > 0:
            print(f"[agent] Failing {len(self.pending_operations)} pending operation(s) for new session")
            for future in self.pending_operations.values():
                if not future.done():
                    future.set_result({
                        "status": "error",
                        "error": "Operation interrupted due to session change. Retry in the new session.",
                    })
            self.pending_operations.clear()

        self.pause_until_user_query = False
        self.executing_cells.clear()
        self.pending_tool_calls.clear()
        self.expected_widgets.clear()
        self.open_stream_blocks.clear()
        self.current_request_id = None

        while not self.pending_messages.empty():
            with suppress(asyncio.QueueEmpty):
                self.pending_messages.get_nowait()

    async def _run_control_query(
        self,
        *,
        prompt: str,
        timeout_submit_seconds: float = 10.0,
        timeout_response_seconds: float = 30.0,
    ) -> ResultMessage | None:
        assert self.client is not None
        session_id = self._current_sdk_session_id()
        await asyncio.wait_for(
            self.client.query(prompt=prompt, session_id=session_id),
            timeout=timeout_submit_seconds,
        )

        receive_iter = self.client.receive_response().__aiter__()
        while True:
            try:
                msg = await asyncio.wait_for(
                    receive_iter.__anext__(),
                    timeout=timeout_response_seconds,
                )
            except StopAsyncIteration:
                return None
            if isinstance(msg, ResultMessage):
                return msg

    async def _handle_stream_event(self, event: dict[str, Any]) -> None:
        event_type = event.get("type")

        if event_type == "message_start":
            self.open_stream_blocks.clear()
            usage = event.get("message", {}).get("usage")
            if isinstance(usage, dict):
                await self._send_usage_update(usage)
            return

        if event_type == "content_block_start":
            block_index = int(event.get("index", -1))
            block = event.get("content_block", {})
            block_type = block.get("type", "unknown")
            if block_index >= 0:
                if block_index in self.open_stream_blocks:
                    await self.send({"type": "agent_stream_block_stop", "block_index": block_index})
                self.open_stream_blocks[block_index] = str(block_type)
            payload: dict[str, object] = {
                "type": "agent_stream_block_start",
                "block_index": block_index,
                "block_type": block_type,
            }
            if block_type == "tool_use":
                payload["block_id"] = str(block.get("id", ""))
                payload["block_name"] = str(block.get("name", ""))
                print(
                    "[agent] tool_use block started "
                    f"index={block_index} name={payload['block_name']} id={payload['block_id']}"
                )
            await self.send(payload)
            return

        if event_type == "content_block_delta":
            block_index = int(event.get("index", -1))
            delta = event.get("delta", {})
            delta_type = delta.get("type")
            block_type_by_delta = {
                "text_delta": "text",
                "thinking_delta": "thinking",
                "input_json_delta": "tool_use",
            }
            if block_index >= 0 and block_index not in self.open_stream_blocks:
                inferred_block_type = block_type_by_delta.get(str(delta_type), "unknown")
                self.open_stream_blocks[block_index] = inferred_block_type
                await self.send({
                    "type": "agent_stream_block_start",
                    "block_index": block_index,
                    "block_type": inferred_block_type,
                })

            if delta_type == "text_delta":
                await self.send({
                    "type": "agent_stream_delta",
                    "block_index": block_index,
                    "block_type": "text",
                    "delta": str(delta.get("text", "")),
                })
            elif delta_type == "thinking_delta":
                await self.send({
                    "type": "agent_stream_delta",
                    "block_index": block_index,
                    "block_type": "thinking",
                    "delta": str(delta.get("thinking", "")),
                })
            elif delta_type == "input_json_delta":
                await self.send({
                    "type": "agent_stream_delta",
                    "block_index": block_index,
                    "block_type": "tool_use",
                    "delta": str(delta.get("partial_json", "")),
                })
            return

        if event_type == "content_block_stop":
            block_index = int(event.get("index", -1))
            if block_index in self.open_stream_blocks:
                self.open_stream_blocks.pop(block_index, None)
            await self.send({
                "type": "agent_stream_block_stop",
                "block_index": block_index,
            })
            return

        if event_type == "message_delta":
            usage = event.get("usage")
            if isinstance(usage, dict):
                await self._send_usage_update(usage)
            return

        if event_type == "message_stop":
            message_usage = event.get("message", {}).get("usage")
            if isinstance(message_usage, dict):
                await self._send_usage_update(message_usage)
            await self._close_open_stream_blocks()
            return

    async def _emit_assistant_fallback(self, msg: AssistantMessage) -> None:
        for idx, block in enumerate(msg.content):
            if isinstance(block, TextBlock):
                await self.send({
                    "type": "agent_stream_block_start",
                    "block_index": idx,
                    "block_type": "text",
                })
                await self.send({
                    "type": "agent_stream_delta",
                    "block_index": idx,
                    "block_type": "text",
                    "delta": block.text,
                })
                await self.send({"type": "agent_stream_block_stop", "block_index": idx})
            elif isinstance(block, ThinkingBlock):
                await self.send({
                    "type": "agent_stream_block_start",
                    "block_index": idx,
                    "block_type": "thinking",
                })
                await self.send({
                    "type": "agent_stream_delta",
                    "block_index": idx,
                    "block_type": "thinking",
                    "delta": block.thinking,
                })
                await self.send({"type": "agent_stream_block_stop", "block_index": idx})
            elif isinstance(block, ToolUseBlock):
                await self.send({
                    "type": "agent_stream_block_start",
                    "block_index": idx,
                    "block_type": "tool_use",
                    "block_id": block.id,
                    "block_name": block.name,
                })
                await self.send({
                    "type": "agent_stream_delta",
                    "block_index": idx,
                    "block_type": "tool_use",
                    "delta": json.dumps(block.input),
                })
                await self.send({"type": "agent_stream_block_stop", "block_index": idx})

    async def _run_query(
        self,
        *,
        prompt: str,
        request_id: str | None,
    ) -> None:
        assert self.client is not None

        self.current_request_id = request_id
        session_id = self._current_sdk_session_id()
        self.open_stream_blocks.clear()

        await self.send({
            "type": "agent_stream_start",
            "timestamp": int(time.time() * 1000),
        })

        try:
            print(
                f"[agent] starting SDK query (request_id={request_id}, session_id={session_id})"
            )
            try:
                await asyncio.wait_for(
                    self.client.query(prompt=prompt, session_id=session_id), timeout=10.0
                )
            except TimeoutError:
                await self.send({
                    "type": "agent_error",
                    "error": "Timed out submitting prompt to Claude runtime",
                    "fatal": False,
                })
                return
            print(f"[agent] SDK query submitted (request_id={request_id})")
            receive_iter = self.client.receive_response().__aiter__()
            while True:
                try:
                    msg = await asyncio.wait_for(receive_iter.__anext__(), timeout=90.0)
                except StopAsyncIteration:
                    break
                except TimeoutError:
                    await self.send({
                        "type": "agent_error",
                        "error": "Timed out waiting for model response from Claude runtime",
                        "fatal": False,
                    })
                    break

                if isinstance(msg, StreamEvent):
                    await self._handle_stream_event(msg.event)
                elif isinstance(msg, AssistantMessage):
                    # Fallback path if partial stream events are disabled.
                    await self._emit_assistant_fallback(msg)
                elif isinstance(msg, ResultMessage):
                    print(
                        "[agent] SDK result message "
                        f"subtype={msg.subtype} is_error={msg.is_error} result={msg.result!r}"
                    )
                    if isinstance(msg.usage, dict):
                        await self._send_usage_update(msg.usage)
                    if msg.is_error:
                        await self.send({
                            "type": "agent_error",
                            "error": msg.result or "Claude query failed",
                            "fatal": False,
                        })
                    # ResultMessage is terminal for the current query.
                    break
                elif isinstance(msg, SystemMessage):
                    continue

            print(f"[agent] finished SDK query (request_id={request_id})")
        except Exception as e:
            await self.send({
                "type": "agent_error",
                "error": f"Agent SDK query failed: {e!s}",
                "fatal": False,
            })
        finally:
            await self._close_open_stream_blocks()
            self.current_request_id = None
            self.current_query_task = None
            await self.send({"type": "agent_stream_complete"})

    async def _run_query_with_turn_prompt(
        self,
        *,
        query: str,
        request_id: str | None,
    ) -> None:
        assert self.client is not None

        try:
            turn_system_prompt = self._compose_turn_system_prompt()
            self.client.options.system_prompt = turn_system_prompt
            turn_prompt = await self._build_turn_prompt(query)
            print(
                "[agent] Turn prompt ready "
                f"(behavior={self.behavior.value}, "
                f"has_plan={self.current_plan is not None}, "
                f"notebook_chars={len(self.latest_notebook_state or '')})"
            )
        except Exception as e:
            await self.send({
                "type": "agent_error",
                "error": f"Failed to build turn prompt: {e!s}",
                "fatal": False,
            })
            self.current_query_task = None
            return

        await self._run_query(prompt=turn_prompt, request_id=request_id)

    # async def _run_quick_inference(self, prompt: str) -> str:
    #     messages = [
    #         {
    #             "role": "user",
    #             "content": prompt
    #             + "\n\nDo not use any markdown formatting. If the content does not contain a clear reasoning process, provide a best-effort summary of the available text. Do not return meta-commentary.",
    #         }
    #     ]

    #     try:
    #         msg = await self.client.messages.create(
    #             model="claude-haiku-4-5-20251001", max_tokens=100, messages=messages
    #         )
    #         return msg.content[0].text

    #     except Exception as e:
    #         print(f"[agent] Failed to run quick inference: {e}")

    #         return ""

    # async def _summarize_and_send_chunk(self, text: str, block_index: int) -> None:
    #     prompt = f"Summarize the most recent thoughts in this reasoning process into a brief, active phrase (2-6 words). Focus on spatial analysis tasks, protocol verification, or scientific reasoning currently being analyzed. If no clear reasoning is present, summarize the general intent. Examples: 'Verifying widget parameters', 'Analyzing QC metrics', 'Checking protocol compliance'.\n\nThinking:\n{text}"
    #     summary = await self._run_quick_inference(prompt)

    #     if summary.strip() != "":
    #         await self.send({
    #             "type": "agent_stream_delta",
    #             "block_index": block_index,
    #             "block_type": "thinking_summary",
    #             "delta": summary,
    #         })

    _context_init_task: asyncio.Task | None = None

    async def handle_init(self, msg: dict[str, object]) -> None:
        new_session_id = self._normalize_session_id(msg.get("session_id"))
        session_changed = new_session_id is not None and new_session_id != self.agent_session_id
        if session_changed:
            print(f"[agent] Session initialized/changed: {new_session_id}")
            self.agent_session_id = new_session_id

        if self.client is not None:
            if session_changed:
                if self.current_query_task is not None and not self.current_query_task.done():
                    print("[agent] Interrupting running query due to session change")
                    await self.client.interrupt()
                    await self._wait_for_running_query_to_stop()
                await self._reset_for_new_session()
            print("[agent] SDK client already initialized; skipping re-init")
            await self.send({"type": "agent_status", "status": "ready"})
            return

        self.init_tools()

        system_prompt_path = context_root.parent / "system_prompt.md"
        self.system_prompt = system_prompt_path.read_text()

        direct_anthropic_key = os.environ.get("AGENT_SDK_DIRECT_ANTHROPIC_KEY", "").strip()
        if direct_anthropic_key != "":
            # Direct mode: bypass Nucleus proxy and call Anthropic directly.
            sdk_env = {
                "ANTHROPIC_AUTH_TOKEN": direct_anthropic_key,
                "ANTHROPIC_API_KEY": direct_anthropic_key,
            }
            print("[agent] SDK gateway mode: direct-anthropic")
        else:
            sdk_base_url = f"{nucleus_url}/infer/plots-agent/anthropic"
            sdk_auth_token = sdk_token if sdk_token != "" else auth_token_sdk
            sdk_env = {
                "ANTHROPIC_BASE_URL": sdk_base_url,
                "ANTHROPIC_AUTH_TOKEN": sdk_auth_token,
                "ANTHROPIC_API_KEY": sdk_auth_token,
            }
            print(f"[agent] SDK gateway mode: nucleus-proxy ({sdk_base_url})")

        def _sdk_stderr(line: str) -> None:
            stripped = line.rstrip()
            if stripped != "":
                print(f"[claude] {stripped}", flush=True)

        self.client = ClaudeSDKClient(
            options=ClaudeAgentOptions(
                system_prompt=self.system_prompt,
                include_partial_messages=True,
                mcp_servers=self.mcp_servers,
                allowed_tools=[
                    "mcp__notebook-tools__create_cell",
                    "mcp__notebook-tools__create_markdown_cell",
                    "mcp__notebook-tools__edit_cell",
                    "mcp__notebook-tools__run_cell",
                    "mcp__notebook-tools__update_plan",
                    "mcp__notebook-tools__submit_response",
                ],
                permission_mode="acceptEdits",
                model="claude-sonnet-4-5",
                max_thinking_tokens=4096,
                env=sdk_env,
                stderr=_sdk_stderr,
            )
        )
        await self.client.connect()
        self.initialized = True

        try:
            mcp_status = await self.client.get_mcp_status()
            print(f"[agent] MCP status: {json.dumps(mcp_status)}")
        except Exception as e:
            print(f"[agent] Failed to fetch MCP status: {e!s}")

        await self.send({"type": "agent_status", "status": "ready"})
        print(
            "[agent] Initialization complete "
            f"(session_id={self._current_sdk_session_id()})"
        )

        # todo(rteqs): planning stuff

    async def handle_query(self, msg: dict[str, object]) -> None:
        assert self.client is not None
        query = msg.get("query", "")
        request_id = msg.get("request_id")
        contextual_node_data = msg.get("contextual_node_data")
        template_version_id = msg.get("template_version_id")
        selected_widgets = msg.get("selected_widgets")
        behavior = msg.get("behavior")

        self.pause_until_user_query = False

        if behavior is not None:
            self.behavior = (
                Behavior.step_by_step
                if behavior == "step_by_step"
                else Behavior.proactive
            )

        full_query = query
        if contextual_node_data:
            full_query = f"{query} \n\nHere is the context of the selected nodes the user would like to use: <ContextualNodeData>{json.dumps(contextual_node_data)}</ContextualNodeData>"

        if self.current_query_task is not None and not self.current_query_task.done():
            await self.send({
                "type": "agent_error",
                "error": "A query is already running. Cancel it before starting another.",
                "fatal": False,
            })
            return

        self.current_query_task = asyncio.create_task(
            self._run_query_with_turn_prompt(query=str(full_query), request_id=request_id)
        )

        print(
            "[agent] Message queued successfully "
            f"(request_id={request_id}, session_id={self._current_sdk_session_id()})"
        )

    async def handle_cancel(self, msg: dict[str, object]) -> None:
        # todo(rteqs): idk if we stil need a request id
        request_id = msg.get("request_id", "unknown")
        print(f"[agent] Cancelling request {request_id}")

        assert self.client is not None

        if self.current_query_task is None or self.current_query_task.done():
            await self.send({
                "type": "agent_error",
                "error": "No active query to cancel.",
                "fatal": False,
            })
            return

        await self.client.interrupt()
        completed_gracefully = await self._wait_for_running_query_to_stop()
        if completed_gracefully:
            print(f"[agent] Cancel completed for request {request_id}")
        else:
            print(f"[agent] Cancel forced task cancellation for request {request_id}")

    async def handle_clear_history(self) -> None:
        assert self.client is not None

        if self.current_query_task is not None and not self.current_query_task.done():
            await self.send({
                "type": "agent_error",
                "error": "Cannot clear history while a query is running. Cancel first.",
                "fatal": False,
            })
            return

        print(f"[agent] Clearing SDK session history (session_id={self._current_sdk_session_id()})")
        try:
            result = await self._run_control_query(prompt="/clear")
            if result is not None:
                print(
                    "[agent] Clear history result "
                    f"subtype={result.subtype} is_error={result.is_error} result={result.result!r}"
                )
                if result.is_error:
                    await self.send({
                        "type": "agent_error",
                        "error": result.result or "Failed to clear history",
                        "fatal": False,
                    })
                    return
        except TimeoutError:
            await self.send({
                "type": "agent_error",
                "error": "Timed out while clearing history.",
                "fatal": False,
            })
            return
        except Exception as e:
            await self.send({
                "type": "agent_error",
                "error": f"Failed to clear history: {e!s}",
                "fatal": False,
            })
            return

        await self.send({
            "type": "agent_history_updated",
            "session_id": self.agent_session_id,
            "request_id": None,
        })

    async def get_full_prompt(self) -> dict:
        self.system_prompt = (context_root.parent / "system_prompt.md").read_text()

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

        # fixme(rteqs): figure out what to do here
        return {
            "system_prompt": self.system_prompt,
            # "messages": messages,
            # "truncated_messages": truncated_messages,
            "model": self.mode_config.get(
                self.mode, ("claude-opus-4-5-20251101", 1024)
            )[0],
            "cells": self.latest_notebook_state
            if self.latest_notebook_state is not None
            else "Interact with agent to populate notebook state.",
            "tree": tree_content,
        }

    async def update_system_prompt(self, msg: dict[str, object]) -> dict:
        assert self.client is not None

        new_content = msg.get("content")
        if not isinstance(new_content, str):
            return {"status": "error", "error": "Invalid content"}

        system_prompt_path = context_root.parent / "system_prompt.md"
        system_prompt_path.write_text(new_content)
        self.system_prompt = new_content
        self.client.options.system_prompt = new_content

        full_prompt = await self.get_full_prompt()
        return {"status": "success", **full_prompt}

    # todo(rteqs): rip a generator for handle_query
    async def accept(self) -> None:
        msg = await self.conn.recv()
        msg_type = msg.get("type")
        msg_request_id = msg.get("request_id")

        print(
            f"[agent] accept: received message type={msg_type} (request_id={msg_request_id})"
        )

        if msg_type == "init":
            print(f"[agent] Message: {msg_type}")
            await self.handle_init(msg)
        elif msg_type == "agent_query":
            query = msg.get("query", "")
            query_preview = query[:60] + "..." if len(query) > 60 else query
            request_id = msg.get("request_id", "unknown")
            print(
                f"[agent] accept: dispatching to handle_query (query={query_preview}, request_id={request_id})"
            )
            handle_start = time.time()
            await self.handle_query(msg)
            handle_elapsed = time.time() - handle_start
            print(
                f"[agent] accept: handle_query completed in {handle_elapsed:.3f}s (request_id={request_id})"
            )
        elif msg_type == "agent_cancel":
            request_id = msg.get("request_id", "unknown")
            print(f"[agent] Cancel: {request_id}")
            await self.handle_cancel(msg)
        elif msg_type == "agent_clear_history":
            print("[agent] Clear history request")
            await self.handle_clear_history()
        elif msg_type == "agent_reset_kernel":
            print("[agent] Reset kernel globals request")
            result = await self.atomic_operation("reset_kernel_globals", {})
            if result.get("status") != "success":
                print(f"[agent] Failed to reset kernel: {result.get('error')}")
        elif msg_type == "agent_action_response":
            print(
                f"[agent] {msg.get('action', 'unknown')} -> {msg.get('status', 'unknown')}"
            )
            await self.handle_action_response(msg)
        elif msg_type == "kernel_message":
            nested_msg = msg.get("message", {})
            nested_type = nested_msg.get("type")
            print(f"[agent] Kernel message: {nested_type}")

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
                if self.pause_until_user_query:
                    print(
                        f"        Suppressing cell {cell_id} result while pause_until_user_query is True"
                    )
                    return
                execution_statuses = {
                    "executing",
                    "awaiting_cell_execution",
                    "thinking",
                    "fixing",
                }
                if (
                    self.current_status is not None
                    and self.current_status not in execution_statuses
                ):
                    print(
                        f"        Not adding cell {cell_id} result because {self.current_status}"
                    )
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
                    print(
                        f"        Cell {cell_id} completed but no active request - updating executing_cells only"
                    )

            elif nested_type == "start_cell":
                cell_id = nested_msg.get("cell_id")
                if cell_id is not None and self.current_request_id is not None:
                    self.executing_cells.add(str(cell_id))
                    print(f"        Added cell {cell_id} to executing_cells")
            elif nested_type == "set_widget_value":
                if self.current_status == "awaiting_user_widget_input":
                    data = nested_msg.get("data", {})
                    for key, value in data.items():
                        if key in self.expected_widgets:
                            self.expected_widgets[key] = value
                            print(f"        Set widget {key}")

                    if all(v is not None for v in self.expected_widgets.values()):
                        self.current_status = "thinking"
                        await self.pending_messages.put({
                            "type": "set_widget_value",
                            "data": self.expected_widgets,
                        })
                        print("        Finished waiting for widget input")
            else:
                print("        Ignored")
        elif msg_type == "get_full_prompt":
            tx_id = msg.get("tx_id")
            print(f"[agent] Get full prompt request (tx_id={tx_id})")

            result = await self.get_full_prompt()
            await self.send({
                "type": "agent_action_response",
                "tx_id": tx_id,
                "status": "success",
                **result,
            })
        elif msg_type == "update_system_prompt":
            tx_id = msg.get("tx_id")
            print(f"[agent] Update system prompt request (tx_id={tx_id})")
            result = await self.update_system_prompt(msg)

            await self.send({"type": "agent_action_response", "tx_id": tx_id, **result})
        elif msg_type == "seed_plan_from_history":
            print("[agent] seed_plan_from_history received")
            plan = msg.get("plan")
            if plan is None or plan.get("steps") is None:
                print("[agent] Invalid plan received")
                return

            if self.current_plan is None:
                self.current_plan = plan
                print(f"[agent] Seeded plan from history: {len(plan['steps'])} steps")
        else:
            print(f"[agent] Unknown message type: {msg_type}")


async def main() -> None:
    global loop
    loop = asyncio.get_running_loop()

    print(f"{datetime.now().isoformat()} [agent] Starting")

    sock = socket.socket(family=socket.AF_UNIX, fileno=int(sys.argv[-1]))
    sock.setblocking(False)

    socket_io_thread = SocketIoThread(socket=sock)
    socket_io_thread.start()
    harness: AgentHarness | None = None
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
        if harness is not None and harness.client is not None:
            try:
                await harness.client.disconnect()
            except Exception:
                traceback.print_exc()
        socket_io_thread.shutdown.set()
        socket_io_thread.join()


if __name__ == "__main__":
    if sys.platform == "linux":
        from ctypes import CDLL

        libc = CDLL("libc.so.6")
        PR_SET_NAME = 15  # https://github.com/torvalds/linux/blob/2df0c02dab829dd89360d98a8a1abaa026ef5798/include/uapi/linux/prctl.h#L56
        libc.prctl(PR_SET_NAME, b"agent")

    asyncio.run(main())
