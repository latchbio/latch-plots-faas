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
    from anthropic import APIStatusError
    from anthropic.types import MessageParam
    from claude_agent_sdk import ClaudeAgentOptions, ClaudeSDKClient
    from claude_agent_sdk.types import (
        AssistantMessage,
        ResultMessage,
        StreamEvent,
        SystemMessage,
        ToolResultBlock,
        ToolUseBlock,
        UserMessage,
    )
except ImportError:
    print("[agent] claude_agent_sdk missing, installing...", flush=True)
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", "claude-agent-sdk"])
    from anthropic import APIStatusError
    from anthropic.types import MessageParam
    from claude_agent_sdk import ClaudeAgentOptions, ClaudeSDKClient
    from claude_agent_sdk.types import (
        AssistantMessage,
        ResultMessage,
        StreamEvent,
        SystemMessage,
        ToolResultBlock,
        ToolUseBlock,
        UserMessage,
    )
from tools import MCP_ALLOWED_TOOL_NAMES, MCP_SERVER_NAME, agent_tools_mcp
from lplots import _inject
from socketio_thread import SocketIoThread
from utils import auth_token_sdk, gql_query, nucleus_url, sdk_token

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
SDK_BUILTIN_ALLOWED_TOOLS = [
    "Read",
    "Grep",
    "Glob",
    "Edit",
    "Write",
    "Bash",
    "WebFetch",
    "WebSearch",
]


class Behavior(Enum):
    proactive = "proactive"
    step_by_step = "step_by_step"


@dataclass
class AgentHarness:
    conn: SocketIoThread
    client: ClaudeSDKClient | None = None
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
    mcp_servers: dict[str, object] = field(
        default_factory=lambda: {MCP_SERVER_NAME: agent_tools_mcp}
    )
    mcp_allowed_tools: list[str] = field(
        default_factory=lambda: list(MCP_ALLOWED_TOOL_NAMES)
    )
    current_query_task: asyncio.Task | None = None
    open_stream_blocks: dict[int, str] = field(default_factory=dict)
    agent_session_metadata: dict[str, object] = field(default_factory=dict)

    async def send(self, msg: dict[str, object]) -> None:
        msg_type = msg.get("type", "unknown")
        if msg_type == "agent_error":
            print(f"[agent] Sending message: agent_error payload={json.dumps(msg)}")
        elif msg_type != "agent_stream_delta":
            print(f"[agent] Sending message: {msg_type}")

        await self.conn.send(msg)

    async def _notify_history_updated(self, *, request_id: str | None = None) -> None:
        payload: dict[str, object] = {
            "type": "agent_history_updated",
            "session_id": self.agent_session_id,
        }
        if request_id is not None:
            payload["request_id"] = request_id
        await self.send(payload)

    async def _fetch_history_from_db(self) -> list[dict[str, object]]:
        if skip_db_history:
            return self.in_memory_history

        if self.agent_session_id is None:
            return []

        resp = await gql_query(
            auth=auth_token_sdk,
            query="""
                query AgentHistory($sessionId: BigInt!) {
                    agentHistories(condition: {sessionId: $sessionId, removed: false}, orderBy: ID_ASC) {
                        nodes { id payload requestId templateVersionId }
                    }
                }
            """,
            variables={"sessionId": str(self.agent_session_id)},
        )

        nodes = resp.get("data", {}).get("agentHistories", {}).get("nodes", [])
        return [{
            "payload": node.get("payload"),
            "request_id": node.get("requestId"),
            "template_version_id": node.get("templateVersionId"),
        } for node in nodes if isinstance(node, dict)]

    async def _build_messages_from_db(self) -> list[MessageParam]:
        history = await self._fetch_history_from_db()
        anthropic_messages: list[MessageParam] = []

        for item in history:
            payload = item.get("payload")
            if not isinstance(payload, dict):
                continue

            payload_type = payload.get("type")
            if payload_type == "anthropic_message":
                role = payload.get("role")
                content = payload.get("content")

                display_widgets = payload.get("display_widgets")
                if role == "user" and display_widgets is not None and isinstance(content, str):
                    for widget in display_widgets:
                        ref_pattern = f"@({widget['widgetKey']}|{widget['id']})"
                        inline_widget = (
                            f"<Widget label=\"{widget['label']}\" type=\"{widget['widgetType']}\" "
                            f"widget_key=\"{widget['widgetKey']}\" cell=\"{widget['cellDisplayName']}\"/>"
                        )
                        content = content.replace(ref_pattern, inline_widget)

                if role == "user" and isinstance(content, dict) and content.get("type") == "cell_result":
                    exception = content.get("exception")
                    logs = content.get("logs")
                    message = content.get("message", "Cell execution completed")

                    content = str(message)
                    if exception:
                        content = f"{content}\n\nException: {exception}"
                    if logs:
                        content = f"{content}\n\nLogs:\n{logs}"

                if isinstance(content, list):
                    cleaned_content: list[object] = []
                    for block in content:
                        if isinstance(block, dict) and block.get("type") == "thinking_summary":
                            continue

                        if isinstance(block, dict) and block.get("type") == "tool_result":
                            block_copy = dict(block)
                            block_content = block_copy.get("content", "{}")
                            if isinstance(block_content, str):
                                try:
                                    result = json.loads(block_content)
                                except Exception:
                                    cleaned_content.append(block_copy)
                                    continue
                                if isinstance(result, dict) and "original_code" in result:
                                    result.pop("original_code", None)
                                    block_copy["content"] = json.dumps(result, sort_keys=True)
                            cleaned_content.append(block_copy)
                            continue

                        cleaned_content.append(block)
                    content = cleaned_content

                template_version_id = item.get("template_version_id")
                if role == "user" and template_version_id is not None:
                    checkpoint_content = (
                        "[auto-generated metadata] "
                        f"template_version_id={template_version_id}"
                    )
                    anthropic_messages.append({"role": "user", "content": checkpoint_content})

                if role in {"user", "assistant"} and isinstance(content, (str, list)):
                    anthropic_messages.append({"role": role, "content": content})
            elif payload_type == "cancellation":
                anthropic_messages.append({
                    "role": "user",
                    "content": "[Request cancelled by user]",
                })

        reordered: list[MessageParam] = []
        pending_tool_ids: set[str] = set()
        deferred: list[MessageParam] = []

        for msg in anthropic_messages:
            role = msg.get("role")
            content = msg.get("content")

            if role == "assistant" and isinstance(content, list):
                for block in content:
                    if isinstance(block, dict) and block.get("type") == "tool_use":
                        tool_id = block.get("id")
                        if tool_id is not None:
                            pending_tool_ids.add(tool_id)

            if len(pending_tool_ids) > 0 and role == "user":
                is_tool_result_message = isinstance(content, list) and any(
                    isinstance(block, dict) and block.get("type") == "tool_result"
                    for block in content
                )
                if is_tool_result_message:
                    reordered.append(msg)
                    if isinstance(content, list):
                        for block in content:
                            if isinstance(block, dict) and block.get("type") == "tool_result":
                                tool_result_id = block.get("tool_use_id")
                                if tool_result_id is not None:
                                    pending_tool_ids.discard(tool_result_id)
                    if len(pending_tool_ids) == 0:
                        reordered.extend(deferred)
                        deferred.clear()
                    continue

                deferred.append(msg)
                continue

            reordered.append(msg)

        reordered.extend(deferred)
        print(f"[agent] Built {len(reordered)} messages from DB")
        return reordered

    async def _insert_history(
        self,
        *,
        event_type: str = "anthropic_message",
        role: str = "user",
        payload: dict[str, object],
        request_id: str | None = None,
        tx_id: str | None = None,
        template_version_id: str | None = None,
    ) -> None:
        event_payload: dict[str, object] = {
            "type": event_type,
            "role": role,
            "timestamp": int(time.time() * 1000),
            **payload,
        }

        if skip_db_history:
            self.in_memory_history.append({
                "payload": event_payload,
                "request_id": request_id,
                "template_version_id": template_version_id,
            })
            await self._notify_history_updated(request_id=request_id)
            return

        if self.agent_session_id is None:
            print("[agent] Skipping history write: missing agent_session_id")
            return

        variables = {
            "sessionId": str(self.agent_session_id),
            "eventType": event_type,
            "payload": event_payload,
            "requestId": request_id,
            "txId": tx_id,
            "templateVersionId": (
                str(template_version_id) if template_version_id is not None else None
            ),
        }

        await gql_query(
            auth=auth_token_sdk,
            query="""
                mutation CreateAgentHistory($sessionId: BigInt!, $eventType: String!, $payload: JSON!, $requestId: String, $txId: String, $templateVersionId: BigInt) {
                    createAgentHistory(input: {agentHistory: {sessionId: $sessionId, eventType: $eventType, payload: $payload, requestId: $requestId, txId: $txId, templateVersionId: $templateVersionId}}) {
                        clientMutationId
                    }
                }
            """,
            variables=variables,
        )

        await self._notify_history_updated(request_id=request_id)

    def _normalize_tool_result_content(
        self, tool_response: object
    ) -> str | list[dict[str, object]]:
        if isinstance(tool_response, str):
            return tool_response

        if isinstance(tool_response, list):
            normalized_blocks: list[dict[str, object]] = []
            for item in tool_response:
                if isinstance(item, dict):
                    block_type = item.get("type")
                    text_value = item.get("text")
                else:
                    block_type = getattr(item, "type", None)
                    text_value = getattr(item, "text", None)
                    if block_type is None and isinstance(text_value, str):
                        block_type = "text"

                if block_type == "text" and isinstance(text_value, str):
                    normalized_blocks.append({
                        "type": "text",
                        "text": text_value,
                    })
                    continue

                return json.dumps(tool_response, default=str)

            if len(normalized_blocks) > 0:
                return normalized_blocks

            return json.dumps(tool_response, default=str)

        if isinstance(tool_response, dict):
            content = tool_response.get("content")
            if isinstance(content, list):
                normalized_content = self._normalize_tool_result_content(content)
                if isinstance(normalized_content, list):
                    return normalized_content
            return json.dumps(tool_response, default=str)

        return json.dumps(tool_response, default=str)

    def _extract_tool_result_text(
        self, content: str | list[dict[str, object]]
    ) -> str | None:
        if isinstance(content, str):
            return content
        for block in content:
            if (
                isinstance(block, dict)
                and block.get("type") == "text"
                and isinstance(block.get("text"), str)
            ):
                return block["text"]
        return None

    def _build_generic_tool_result_content(
        self,
        *,
        tool_name: str,
        content: str | list[dict[str, object]],
        is_error: bool,
    ) -> str | list[dict[str, object]]:
        text_content = self._extract_tool_result_text(content)
        if text_content is None:
            return content

        try:
            parsed = json.loads(text_content)
        except Exception:
            parsed = None
        if isinstance(parsed, dict) and isinstance(parsed.get("tool_name"), str):
            return content

        result_payload: dict[str, object] = {
            "tool_name": tool_name,
            "success": not is_error,
            "summary": (
                f"{tool_name} failed"
                if is_error
                else f"{tool_name} completed"
            ),
        }
        if text_content.strip() != "":
            result_payload["raw_output"] = text_content
        if is_error:
            result_payload["error"] = text_content
            if text_content.strip() != "":
                result_payload["summary"] = f"{tool_name} failed: {text_content[:200]}"

        return [{"type": "text", "text": json.dumps(result_payload)}]

    def _extract_tool_blocks_from_sdk_message(
        self, msg: object
    ) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
        content = getattr(msg, "content", None)
        if not isinstance(content, list):
            return [], []

        tool_use_blocks: list[dict[str, object]] = []
        tool_result_blocks: list[dict[str, object]] = []

        for block in content:
            if isinstance(block, ToolUseBlock):
                tool_use_id = str(block.id).strip()
                tool_name = str(block.name).strip()
                if tool_use_id == "" or tool_name == "":
                    continue
                tool_input = block.input if isinstance(block.input, dict) else {}
                tool_use_blocks.append({
                    "type": "tool_use",
                    "id": tool_use_id,
                    "name": tool_name,
                    "input": tool_input,
                })
                continue

            if not isinstance(block, ToolResultBlock):
                continue

            tool_use_id = str(block.tool_use_id).strip()
            if tool_use_id == "":
                continue
            normalized_content = self._normalize_tool_result_content(block.content)
            payload_block: dict[str, object] = {
                "type": "tool_result",
                "tool_use_id": tool_use_id,
                "content": normalized_content,
            }
            if isinstance(block.is_error, bool):
                payload_block["is_error"] = block.is_error
            tool_result_blocks.append(payload_block)

        return tool_use_blocks, tool_result_blocks

    async def _persist_tool_blocks_from_sdk_message(
        self,
        *,
        msg: object,
        request_id: str | None,
        tool_use_index: dict[str, str],
    ) -> None:
        tool_use_blocks, tool_result_blocks = self._extract_tool_blocks_from_sdk_message(msg)
        message_role: str | None = None
        if isinstance(msg, AssistantMessage):
            message_role = "assistant"
        elif isinstance(msg, UserMessage):
            message_role = "user"

        for block in tool_use_blocks:
            tool_use_id = block.get("id")
            tool_name = block.get("name")
            if not isinstance(tool_use_id, str) or tool_use_id == "":
                continue
            if not isinstance(tool_name, str) or tool_name == "":
                continue
            tool_use_index[tool_use_id] = tool_name

        if len(tool_use_blocks) > 0:
            try:
                await self._insert_history(
                    role=message_role or "assistant",
                    payload={"content": tool_use_blocks},
                    request_id=request_id,
                )
            except Exception as e:
                print(f"[agent] Failed to persist message tool_use blocks: {e!s}")

        if len(tool_result_blocks) > 0:
            tool_result_blocks_with_added_fields: list[dict[str, object]] = []
            for block in tool_result_blocks:
                tool_result_block_with_added_fields = dict(block)
                tool_use_id = block.get("tool_use_id")
                if isinstance(tool_use_id, str):
                    tool_name = tool_use_index.get(tool_use_id)
                    content = block.get("content")
                    if (
                        isinstance(tool_name, str)
                        and (
                            isinstance(content, str)
                            or isinstance(content, list)
                        )
                    ):
                        tool_result_block_with_added_fields["content"] = (
                            self._build_generic_tool_result_content(
                                tool_name=tool_name,
                                content=content,
                                is_error=bool(block.get("is_error", False)),
                            )
                        )
                tool_result_blocks_with_added_fields.append(
                    tool_result_block_with_added_fields
                )

            try:
                await self._insert_history(
                    role=message_role or "user",
                    payload={"content": tool_result_blocks_with_added_fields},
                    request_id=request_id,
                )
            except Exception as e:
                print(f"[agent] Failed to persist message tool_result blocks: {e!s}")

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
        behavior_file = self._behavior_file_name()
        turn_behavior_content = (context_root / "turn_behavior" / behavior_file).read_text()
        examples_content = (context_root / "examples" / behavior_file).read_text()

        notebook_state = await self.refresh_cells_context()
        self.latest_notebook_state = notebook_state

        context_blocks = [
            f"<turn_behavior>\n{turn_behavior_content}\n</turn_behavior>",
            f"<examples>\n{examples_content}\n</examples>",
            f"<current_notebook_state>\n{notebook_state}\n</current_notebook_state>",
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

    def _normalize_session_id(self, raw_session_id: object) -> str | None:
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

    async def _stop_conversation_loop(self) -> None:
        self.conversation_running = False
        if self.conversation_task is None:
            return
        if not self.conversation_task.done():
            self.conversation_task.cancel()
            with suppress(asyncio.CancelledError):
                await self.conversation_task
        self.conversation_task = None

    def _start_conversation_loop(self) -> None:
        if self.conversation_task is not None and not self.conversation_task.done():
            return
        self.conversation_task = asyncio.create_task(self.run_agent_loop())

    async def _wait_for_message(self) -> tuple[str, str | None]:
        print("[agent] _wait_for_message: waiting for message...")
        while True:
            msg = await self.pending_messages.get()
            msg_type = msg.get("type")

            if (
                self.pause_until_user_query
                and msg_type in {"resume", "cell_result", "set_widget_value"}
            ):
                print(f"[agent] Suppressing follow-up message type={msg_type}")
                continue

            if msg_type == "resume":
                request_id = msg.get("request_id")
                if isinstance(request_id, str):
                    self.current_request_id = request_id
                self.current_status = "thinking"
                content = msg.get("content")
                if isinstance(content, str) and content.strip() != "":
                    return content, self.current_request_id
                return "Continue with the next step.", self.current_request_id

            if msg_type == "user_query":
                request_id = msg.get("request_id")
                self.current_request_id = request_id if isinstance(request_id, str) else None
                self.current_status = "thinking"
                self.pause_until_user_query = False

                content = msg.get("content")
                if isinstance(content, str) and content.strip() != "":
                    return content, self.current_request_id
                print("[agent] Ignoring user_query with empty content")
                continue

            if msg_type == "cell_result":
                if self.current_request_id is None:
                    print(
                        f"[agent] Ignoring cell_result for {msg.get('cell_id')} - no active request"
                    )
                    continue

                cell_id = msg.get("cell_id")
                success = bool(msg.get("success", True))
                cell_name = msg.get("display_name")
                logs = msg.get("logs")

                if success:
                    result_message = f"Cell {cell_name} ({cell_id}) executed successfully."
                else:
                    exception = msg.get("exception", "Unknown error")
                    result_message = (
                        f"Cell {cell_name} ({cell_id}) execution failed. Exception: {exception}"
                    )

                try:
                    await self._insert_history(
                        role="user",
                        payload={
                            "content": result_message,
                            "hidden": True,
                        },
                        request_id=self.current_request_id,
                    )
                except Exception as e:
                    print(f"[agent] Failed to persist cell_result history: {e!s}")

                logs_text = ""
                if isinstance(logs, str) and logs.strip() != "":
                    logs_text = f"\nRecent logs:\n{logs}"

                if self.pending_auto_continue and len(self.executing_cells) == 0:
                    self.pending_auto_continue = False
                    return (
                        "A requested cell finished running. Continue with the next step.\n"
                        + result_message
                        + logs_text,
                        self.current_request_id,
                    )

                return (
                    "Cell execution update:\n" + result_message + logs_text,
                    self.current_request_id,
                )

            if msg_type == "set_widget_value":
                if self.current_request_id is None:
                    print("[agent] Ignoring set_widget_value - no active request")
                    continue

                data = msg.get("data", {})
                if not isinstance(data, dict):
                    data = {}
                widget_info = ", ".join(f"{k}={v}" for k, v in data.items())
                content = f"User provided input via widget(s): {widget_info}"

                try:
                    await self._insert_history(
                        role="user",
                        payload={
                            "content": content,
                            "hidden": True,
                        },
                        request_id=self.current_request_id,
                    )
                except Exception as e:
                    print(f"[agent] Failed to persist widget history: {e!s}")

                return content, self.current_request_id

            print(f"[agent] Unknown pending message type={msg_type}, ignoring")

    async def _complete_turn(self) -> None:
        if self.current_request_id is None:
            return

        should_continue = self.should_auto_continue
        self.should_auto_continue = False

        await self._notify_history_updated(request_id=self.current_request_id)

        if should_continue:
            print("[agent] Auto-continuing as requested by model")
            await self.pending_messages.put({
                "type": "resume",
                "content": "Continue with the next step.",
                "request_id": self.current_request_id,
            })

    async def run_agent_loop(self) -> None:
        print("[agent] run_agent_loop: started")
        self.conversation_running = True

        while self.conversation_running:
            try:
                query, request_id = await self._wait_for_message()
            except asyncio.CancelledError:
                break
            except Exception as e:
                print(f"[agent] _wait_for_message failed: {e!s}")
                continue

            if not self.conversation_running:
                break

            turn_task: asyncio.Task | None = None
            try:
                turn_task = asyncio.create_task(
                    self._run_query_with_turn_prompt(query=query, request_id=request_id)
                )
                self.current_query_task = turn_task
                await turn_task
            except Exception as e:
                print(f"[agent] run_agent_loop turn failed: {e!s}")
            finally:
                if turn_task is not None and self.current_query_task is turn_task:
                    self.current_query_task = None
                await self._complete_turn()

        print("[agent] run_agent_loop: stopped")

    async def _disconnect_sdk_client(self) -> None:
        if self.client is None:
            return
        try:
            await self.client.disconnect()
        finally:
            self.client = None

    async def _load_agent_session_metadata(self, session_id: str) -> dict[str, object]:
        response = await gql_query(
            auth=auth_token_sdk,
            query="""
                query AgentSessionMetadata($id: BigInt!) {
                    agentSession(id: $id) {
                        metadata
                    }
                }
            """,
            variables={"id": session_id},
        )

        data = response.get("data")
        if not isinstance(data, dict):
            return {}
        session = data.get("agentSession")
        if not isinstance(session, dict):
            return {}
        metadata = session.get("metadata")
        if not isinstance(metadata, dict):
            return {}
        return metadata

    async def _persist_agent_session_metadata(self, metadata: dict[str, object]) -> None:
        if self.agent_session_id is None:
            return

        await gql_query(
            auth=auth_token_sdk,
            query="""
                mutation UpdateAgentSessionMetadata($id: BigInt!, $metadata: JSON) {
                    updateAgentSession(input: {id: $id, patch: {metadata: $metadata}}) {
                        clientMutationId
                    }
                }
            """,
            variables={
                "id": self.agent_session_id,
                "metadata": metadata,
            },
        )

    async def _capture_claude_session_id(self, raw_session_id: object) -> None:
        session_id = self._normalize_session_id(raw_session_id)
        if session_id is None or self.agent_session_id is None:
            return

        existing = self._normalize_session_id(
            self.agent_session_metadata.get("claude_session_id")
        )
        if existing == session_id:
            return

        merged_metadata = dict(self.agent_session_metadata)
        merged_metadata["claude_session_id"] = session_id

        try:
            await self._persist_agent_session_metadata(merged_metadata)
        except Exception as e:
            print(f"[agent] Failed to persist Claude session id metadata: {e!s}")
            return

        self.agent_session_metadata = merged_metadata
        print(
            "[agent] Persisted Claude session id "
            f"(db_session_id={self.agent_session_id}, claude_session_id={session_id})"
        )

    def _build_sdk_env(self) -> dict[str, str]:
        direct_anthropic_key = os.environ.get("AGENT_SDK_DIRECT_ANTHROPIC_KEY", "").strip()
        if direct_anthropic_key != "":
            print("[agent] SDK gateway mode: direct-anthropic")
            return {
                "ANTHROPIC_AUTH_TOKEN": direct_anthropic_key,
                "ANTHROPIC_API_KEY": direct_anthropic_key,
            }

        sdk_base_url = f"{nucleus_url}/infer/plots-agent/anthropic"
        sdk_auth_token = sdk_token if sdk_token != "" else auth_token_sdk
        print(f"[agent] SDK gateway mode: nucleus-proxy ({sdk_base_url})")
        return {
            "ANTHROPIC_BASE_URL": sdk_base_url,
            "ANTHROPIC_AUTH_TOKEN": sdk_auth_token,
            "ANTHROPIC_API_KEY": sdk_auth_token,
        }

    async def _connect_sdk_client(self, *, resume_session_id: str | None) -> None:
        self.system_prompt = self._compose_turn_system_prompt()
        sdk_env = self._build_sdk_env()

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
                    *self.mcp_allowed_tools,
                    *SDK_BUILTIN_ALLOWED_TOOLS,
                ],
                permission_mode="acceptEdits",
                model="claude-opus-4-6",
                thinking={"type": "adaptive"},
                resume=resume_session_id,
                env=sdk_env,
                stderr=_sdk_stderr,
            )
        )
        await self.client.connect()

        try:
            mcp_status = await self.client.get_mcp_status()
            print(f"[agent] MCP status: {json.dumps(mcp_status)}")
        except Exception as e:
            print(f"[agent] Failed to fetch MCP status: {e!s}")

    async def _reset_for_new_session(self) -> None:
        await self._stop_conversation_loop()

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
        self.current_status = None
        self.should_auto_continue = False
        self.pending_auto_continue = False

        while not self.pending_messages.empty():
            with suppress(asyncio.QueueEmpty):
                self.pending_messages.get_nowait()

        self.agent_session_metadata = {}

    async def _handle_stream_event(
        self,
        event: dict[str, Any],
        *,
        collected_assistant_blocks: dict[int, dict[str, object]] | None = None,
    ) -> None:
        event_type = event.get("type")

        if event_type == "message_start":
            self.open_stream_blocks.clear()
            return

        if event_type == "content_block_start":
            block_index = int(event.get("index", -1))
            block = event.get("content_block", {})
            block_type = block.get("type", "unknown")
            if block_index >= 0:
                if block_index in self.open_stream_blocks:
                    await self.send({"type": "agent_stream_block_stop", "block_index": block_index})
                self.open_stream_blocks[block_index] = str(block_type)

            if collected_assistant_blocks is not None and block_index >= 0:
                if block_type == "text":
                    collected_assistant_blocks[block_index] = {"type": "text", "text": ""}
                elif block_type == "thinking":
                    collected_assistant_blocks[block_index] = {
                        "type": "thinking",
                        "thinking": "",
                    }

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

            if collected_assistant_blocks is not None and block_index >= 0:
                existing = collected_assistant_blocks.get(block_index)
                if existing is not None:
                    if delta_type == "text_delta":
                        existing["text"] = str(existing.get("text", "")) + str(delta.get("text", ""))
                    elif delta_type == "thinking_delta":
                        existing["thinking"] = str(existing.get("thinking", "")) + str(
                            delta.get("thinking", "")
                        )
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
            return

        if event_type == "message_stop":
            await self._close_open_stream_blocks()
            return

    async def _run_query(
        self,
        *,
        prompt: str,
        request_id: str | None,
    ) -> None:
        assert self.client is not None

        session_id = self._current_sdk_session_id()
        self.open_stream_blocks.clear()
        run_started_at = time.perf_counter()
        assistant_blocks_by_index: dict[int, dict[str, object]] = {}
        assistant_message_started_at: float | None = None
        persisted_assistant_message_this_turn = False
        terminal_error: str | None = None
        stream_complete_error: dict[str, object] | None = None
        usage_data: dict[str, Any] | None = None
        tool_use_index: dict[str, str] = {}
        assistant_error_type: str | None = None

        async def persist_current_assistant_message(
            *, duration_seconds: float | None = None
        ) -> bool:
            nonlocal persisted_assistant_message_this_turn
            if len(assistant_blocks_by_index) == 0:
                return False

            assistant_content: list[dict[str, object]] = []
            for idx in sorted(assistant_blocks_by_index):
                assistant_content.append(assistant_blocks_by_index[idx])

            resolved_duration = duration_seconds
            if resolved_duration is None:
                duration_base = (
                    assistant_message_started_at
                    if assistant_message_started_at is not None
                    else run_started_at
                )
                resolved_duration = max(0.0, time.perf_counter() - duration_base)

            try:
                await self._insert_history(
                    role="assistant",
                    payload={
                        "content": assistant_content,
                        "duration": resolved_duration,
                    },
                    request_id=request_id,
                )
                persisted_assistant_message_this_turn = True
            except Exception as e:
                print(f"[agent] Failed to persist assistant history: {e!s}")
            finally:
                assistant_blocks_by_index.clear()

            return True

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
                terminal_error = "Timed out submitting prompt to Claude runtime"
                await self.send({
                    "type": "agent_error",
                    "error": terminal_error,
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
                    terminal_error = "Timed out waiting for model response from Claude runtime"
                    break

                if isinstance(msg, StreamEvent):
                    await self._capture_claude_session_id(msg.session_id)
                    event_type = msg.event.get("type")
                    if event_type == "message_start":
                        assistant_blocks_by_index.clear()
                        assistant_message_started_at = time.perf_counter()
                        message_usage = msg.event.get("message", {}).get("usage")
                        if isinstance(message_usage, dict):
                            usage_data = message_usage
                    elif event_type == "message_delta":
                        delta_usage = msg.event.get("usage")
                        if isinstance(delta_usage, dict):
                            usage_data = delta_usage
                    elif event_type == "message_stop":
                        stop_usage = msg.event.get("message", {}).get("usage")
                        if isinstance(stop_usage, dict):
                            usage_data = stop_usage
                    await self._handle_stream_event(
                        msg.event, collected_assistant_blocks=assistant_blocks_by_index
                    )
                    if event_type == "message_stop":
                        message_duration_seconds: float | None = None
                        if assistant_message_started_at is not None:
                            message_duration_seconds = max(
                                0.0, time.perf_counter() - assistant_message_started_at
                            )
                        await persist_current_assistant_message(
                            duration_seconds=message_duration_seconds
                        )
                        assistant_message_started_at = None
                elif isinstance(msg, ResultMessage):
                    await self._capture_claude_session_id(msg.session_id)
                    print(
                        "[agent] SDK result message "
                        f"subtype={msg.subtype} is_error={msg.is_error} result={msg.result!r}"
                    )
                    if msg.is_error:
                        terminal_error = msg.result if msg.result is not None else "Claude query failed"
                        if (
                            isinstance(msg.result, str)
                            and "prompt is too long" in msg.result.lower()
                        ):
                            stream_complete_error = {
                                "message": (
                                    "This conversation is too long for the agent to continue. "
                                    "Clear history and try again."
                                ),
                                "should_contact_support": False,
                                "should_clear_history": True,
                            }
                        elif assistant_error_type in {
                            "authentication_failed",
                            "billing_error",
                            "invalid_request",
                        }:
                            stream_complete_error = {
                                "message": "An unexpected error occurred. Please try again.",
                                "should_contact_support": True,
                                "should_clear_history": False,
                            }
                        elif assistant_error_type in {
                            "rate_limit",
                            "server_error",
                            "unknown",
                        }:
                            stream_complete_error = {
                                "message": (
                                    "Our model provider is experiencing a temporary issue. "
                                    "Please try again in a few minutes."
                                ),
                                "should_contact_support": False,
                                "should_clear_history": False,
                            }

                        if stream_complete_error is not None:
                            agent_error_message = str(stream_complete_error["message"])
                        else:
                            agent_error_message = terminal_error
                        await self.send({
                            "type": "agent_error",
                            "error": agent_error_message,
                            "fatal": False,
                        })
                    # todo(tim): reconsider if needed 
                    elif (
                        len(assistant_blocks_by_index) == 0
                        and not persisted_assistant_message_this_turn
                        and isinstance(msg.result, str)
                        and msg.result != ""
                    ):
                        assistant_blocks_by_index[0] = {
                            "type": "text",
                            "text": msg.result,
                        }
                        if assistant_message_started_at is None:
                            assistant_message_started_at = run_started_at
                    break
                elif isinstance(msg, SystemMessage):
                    if msg.subtype == "init":
                        await self._capture_claude_session_id(msg.data.get("session_id"))
                    continue
                else:
                    if isinstance(msg, AssistantMessage) and msg.error is not None:
                        assistant_error_type = str(msg.error)
                        print(f"[agent] assistant message error={msg.error}")
                    await self._capture_claude_session_id(getattr(msg, "session_id", None))
                    await self._persist_tool_blocks_from_sdk_message(
                        msg=msg, request_id=request_id, tool_use_index=tool_use_index
                    )
                    continue

            print(f"[agent] finished SDK query (request_id={request_id})")
        except APIStatusError as e:
            should_contact_support = 400 <= e.status_code < 500 and e.status_code != 429

            if should_contact_support:
                user_message = "An unexpected error occurred. Please try again."
            else:
                user_message = (
                    "Our model provider is experiencing a temporary issue. "
                    "Please try again in a few minutes."
                )

            stream_complete_error = {
                "message": user_message,
                "should_contact_support": should_contact_support,
                "should_clear_history": False,
            }
            terminal_error = f"API error: {e!s}"
            await self.send({
                "type": "agent_error",
                "error": terminal_error,
                "fatal": False,
            })
        except Exception as e:
            terminal_error = f"Agent SDK query failed: {e!s}"
            await self.send({
                "type": "agent_error",
                "error": terminal_error,
                "fatal": False,
            })
        finally:
            pending_duration_base = (
                assistant_message_started_at
                if assistant_message_started_at is not None
                else run_started_at
            )
            await persist_current_assistant_message(
                duration_seconds=max(0.0, time.perf_counter() - pending_duration_base)
            )

            if terminal_error is not None:
                try:
                    if stream_complete_error is not None:
                        error_history_payload = dict(stream_complete_error)
                        error_history_payload["raw_error"] = terminal_error
                        if assistant_error_type is not None:
                            error_history_payload["assistant_error_type"] = assistant_error_type
                    else:
                        error_history_payload = {
                            "message": terminal_error,
                            "should_contact_support": False,
                        }
                    await self._insert_history(
                        event_type="error",
                        role="system",
                        payload=error_history_payload,
                        request_id=request_id,
                    )
                except Exception as e:
                    print(f"[agent] Failed to persist error history: {e!s}")

            await self._close_open_stream_blocks()
            if usage_data is not None:
                await self._send_usage_update(usage_data)
            stream_complete_payload: dict[str, object] = {"type": "agent_stream_complete"}
            if stream_complete_error is not None:
                stream_complete_payload["error"] = stream_complete_error
            await self.send(stream_complete_payload)

    async def _run_query_with_turn_prompt(
        self,
        *,
        query: str,
        request_id: str | None,
    ) -> None:
        assert self.client is not None

        try:
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
        elif self.agent_session_id is None:
            self.agent_session_id = new_session_id

        if self.agent_session_id is None:
            print("[agent] Missing session_id in init; cannot initialize Claude SDK client")
            await self.send({
                "type": "agent_error",
                "error": "Missing agent session id.",
                "fatal": False,
            })
            return

        if self.client is not None and not session_changed:
            print("[agent] SDK client already initialized; skipping re-init")
            await self.send({"type": "agent_status", "status": "ready"})
            self._start_conversation_loop()
            return

        if self.client is not None and session_changed:
            if self.current_query_task is not None and not self.current_query_task.done():
                print("[agent] Interrupting running query due to session change")
                await self.client.interrupt()
                await self._wait_for_running_query_to_stop()
            await self._reset_for_new_session()
            await self._disconnect_sdk_client()

        try:
            metadata = await self._load_agent_session_metadata(self.agent_session_id)
        except Exception as e:
            print(f"[agent] Failed to load agent session metadata: {e!s}")
            metadata = {}
        self.agent_session_metadata = metadata
        resume_session_id = self._normalize_session_id(metadata.get("claude_session_id"))

        if self.client is None:
            await self._connect_sdk_client(resume_session_id=resume_session_id)
            if resume_session_id is None:
                print("[agent] SDK initialized without resume session")
            else:
                print(f"[agent] SDK initialized with resume session {resume_session_id}")

        await self.send({"type": "agent_status", "status": "ready"})
        print(
            "[agent] Initialization complete "
            f"(db_session_id={self.agent_session_id}, "
            f"claude_session_id={resume_session_id})"
        )
        self._start_conversation_loop()

        # todo(rteqs): planning stuff

    async def handle_query(self, msg: dict[str, object]) -> None:
        assert self.client is not None
        query = msg.get("query", "")
        raw_request_id = msg.get("request_id")
        request_id = str(raw_request_id) if raw_request_id is not None else None
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

        try:
            user_payload: dict[str, object] = {
                "content": str(full_query),
            }
            if query is not None:
                user_payload["display_query"] = query
            if contextual_node_data is not None:
                user_payload["display_nodes"] = contextual_node_data
            if selected_widgets is not None:
                user_payload["display_widgets"] = selected_widgets
            await self._insert_history(
                role="user",
                payload=user_payload,
                request_id=request_id if isinstance(request_id, str) else None,
                template_version_id=(
                    str(template_version_id)
                    if template_version_id is not None
                    else None
                ),
            )
        except Exception as e:
            print(f"[agent] Failed to persist user history: {e!s}")

        await self.pending_messages.put({
            "type": "user_query",
            "content": str(full_query),
            "request_id": request_id if isinstance(request_id, str) else None,
        })

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

        try:
            await self._insert_history(
                event_type="cancellation",
                role="system",
                payload={"reason": "Request cancelled by user"},
                request_id=str(request_id),
            )
        except Exception as e:
            print(f"[agent] Failed to persist cancellation history: {e!s}")

        self.should_auto_continue = False
        self.pending_auto_continue = False
        self.pause_until_user_query = True
        self.expected_widgets.clear()
        self.current_status = None
        self.executing_cells.clear()

    async def get_full_prompt(self) -> dict:
        self.system_prompt = (context_root.parent / "system_prompt.md").read_text()
        messages = await self._build_messages_from_db()

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

        return {
            "system_prompt": self.system_prompt,
            "messages": messages,
            "model": "claude-opus-4-6",
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
        if harness is not None:
            try:
                await harness._disconnect_sdk_client()
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