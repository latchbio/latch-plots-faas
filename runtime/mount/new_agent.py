import asyncio
import json
import os
import socket
import sys
import time
import traceback
import uuid
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any

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
from utils import auth_token_sdk, nucleus_url

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
    agent_session_id: int | None = None
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

    mode_config: dict[Mode, tuple[str, int | None]] = field(
        default_factory=lambda: {
            Mode.planning: ("claude-opus-4-5-20251101", 4096),
            Mode.executing: ("claude-opus-4-5-20251101", 1024),
            Mode.debugging: ("claude-opus-4-5-20251101", 2048),
        }
    )

    async def send(self, msg: dict[str, object]) -> None:
        msg_type = msg.get("type", "unknown")
        if msg_type != "agent_stream_delta":
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

    async def _handle_stream_event(self, event: dict[str, Any]) -> None:
        event_type = event.get("type")

        if event_type == "message_start":
            usage = event.get("message", {}).get("usage")
            if isinstance(usage, dict):
                await self._send_usage_update(usage)
            return

        if event_type == "content_block_start":
            block_index = int(event.get("index", -1))
            block = event.get("content_block", {})
            block_type = block.get("type", "unknown")
            payload: dict[str, object] = {
                "type": "agent_stream_block_start",
                "block_index": block_index,
                "block_type": block_type,
            }
            if block_type == "tool_use":
                payload["block_id"] = str(block.get("id", ""))
                payload["block_name"] = str(block.get("name", ""))
            await self.send(payload)
            return

        if event_type == "content_block_delta":
            block_index = int(event.get("index", -1))
            delta = event.get("delta", {})
            delta_type = delta.get("type")
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
            await self.send({
                "type": "agent_stream_block_stop",
                "block_index": int(event.get("index", -1)),
            })
            return

        if event_type == "message_delta":
            usage = event.get("usage")
            if isinstance(usage, dict):
                await self._send_usage_update(usage)
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
        session_id = (
            str(self.agent_session_id) if self.agent_session_id is not None else "default"
        )

        await self.send({
            "type": "agent_stream_start",
            "timestamp": int(time.time() * 1000),
        })

        try:
            await self.client.query(prompt=prompt, session_id=session_id)
            async for msg in self.client.receive_response():
                if isinstance(msg, StreamEvent):
                    await self._handle_stream_event(msg.event)
                elif isinstance(msg, AssistantMessage):
                    # Fallback path if partial stream events are disabled.
                    await self._emit_assistant_fallback(msg)
                elif isinstance(msg, ResultMessage):
                    if isinstance(msg.usage, dict):
                        await self._send_usage_update(msg.usage)
                    if msg.is_error:
                        await self.send({
                            "type": "agent_error",
                            "error": msg.result or "Claude query failed",
                            "fatal": False,
                        })
                elif isinstance(msg, SystemMessage):
                    continue
        except Exception as e:
            await self.send({
                "type": "agent_error",
                "error": f"Agent SDK query failed: {e!s}",
                "fatal": False,
            })
        finally:
            self.current_request_id = None
            self.current_query_task = None
            await self.send({"type": "agent_stream_complete"})

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
        self.init_tools()

        system_prompt_path = context_root.parent / "system_prompt.md"
        self.system_prompt = system_prompt_path.read_text()

        sdk_base_url = f"{nucleus_url}/infer/plots-agent/anthropic"
        sdk_auth_token = auth_token_sdk
        sdk_env = {
            "ANTHROPIC_BASE_URL": sdk_base_url,
            "ANTHROPIC_AUTH_TOKEN": sdk_auth_token,
            "ANTHROPIC_API_KEY": sdk_auth_token,
        }

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
                ],
                permission_mode="bypassPermissions",
                model="claude-opus-4-5-20251101",
                max_thinking_tokens=4096,
                env=sdk_env,
            )
        )
        await self.client.connect()

        await self.send({"type": "agent_status", "status": "ready"})
        print("[agent] Initialization complete")

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
            self._run_query(prompt=str(full_query), request_id=request_id)
        )

        print(f"[agent] Message queued successfully (request_id={request_id})")

    async def handle_cancel(self, msg: dict[str, object]) -> None:
        # todo(rteqs): idk if we stil need a request id
        request_id = msg.get("request_id", "unknown")
        print(f"[agent] Cancelling request {request_id}")

        assert self.client is not None

        await self.client.interrupt()

        if self.current_query_task is not None and not self.current_query_task.done():
            try:
                await asyncio.wait_for(self.current_query_task, timeout=2.0)
            except TimeoutError:
                self.current_query_task.cancel()
            except asyncio.CancelledError:
                pass

    async def handle_clear_history(self) -> None:
        assert self.client is not None

        await self.client.query(prompt="/clear")
        await self.send({
            "type": "agent_history_updated",
            "session_id": str(self.agent_session_id)
            if self.agent_session_id is not None
            else None,
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
