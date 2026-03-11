import asyncio
import json
import os
import socket
import sys
import time
import traceback
import uuid
from collections.abc import Iterable
from contextlib import suppress
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Literal, TypedDict

from anthropic.types import (
    MessageDeltaUsage,
    MessageParam,
    RawContentBlockDeltaEvent,
    RawContentBlockStartEvent,
    RawContentBlockStopEvent,
    RawMessageDeltaEvent,
    RawMessageStartEvent,
    RawMessageStopEvent,
    Usage,
)
from claude_agent_sdk import ClaudeAgentOptions, ClaudeSDKClient, HookMatcher
from claude_agent_sdk.types import (
    AssistantMessage,
    HookContext,
    McpSdkServerConfig,
    PostToolUseHookInput,
    ResultMessage,
    StreamEvent,
    SyncHookJSONOutput,
    SystemMessage,
    ToolResultBlock,
    ToolUseBlock,
    UserMessage,
)
from latch_data_validation.data_validation import validate
from lplots import _inject
from socketio_thread import SocketIoThread
from tools import MCP_ALLOWED_TOOL_NAMES, MCP_SERVER_NAME, agent_tools_mcp
from utils import auth_token_sdk, gql_query, nucleus_url, pod_id, sdk_token

sys.stdout.reconfigure(line_buffering=True)

skip_db_history = os.environ.get("AGENT_SKIP_DB_HISTORY") == "1"
legacy_history_window_size = 20

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
NOTEBOOK_MUTATION_TOOL_NAMES = (
    "create_cell",
    "create_markdown_cell",
    "edit_cell",
    "delete_cell",
    "delete_all_cells",
    "rename_notebook",
    "create_tab",
    "rename_tab",
    "restore_checkpoint",
)
NOTEBOOK_MUTATION_TOOL_MATCHER = (
    "^mcp__" + MCP_SERVER_NAME + "__(" + "|".join(NOTEBOOK_MUTATION_TOOL_NAMES) + ")$"
)


class Behavior(Enum):
    proactive = "proactive"
    step_by_step = "step_by_step"


HistoryRole = Literal["user", "assistant", "system"]


class AgentSessionQueryNode(TypedDict):
    claudeSessionId: str | None


class AgentSessionQueryData(TypedDict):
    agentSession: AgentSessionQueryNode | None


class AgentSessionQueryResp(TypedDict):
    data: AgentSessionQueryData


class TextToolResultBlock(TypedDict):
    type: Literal["text"]
    text: str


class ImageToolResultSource(TypedDict):
    type: Literal["base64"]
    media_type: str
    data: str


class ImageToolResultBlock(TypedDict):
    type: Literal["image"]
    source: ImageToolResultSource


ToolResultContent = str | list[TextToolResultBlock | ImageToolResultBlock]


class ToolUseHistoryBlock(TypedDict):
    type: Literal["tool_use"]
    id: str
    name: str
    input: dict[str, Any]


class ToolResultHistoryBlock(TypedDict, total=False):
    type: Literal["tool_result"]
    tool_use_id: str
    content: ToolResultContent
    is_error: bool


class GenericToolResultPayload(TypedDict, total=False):
    tool_name: str
    success: bool
    summary: str
    raw_output: str
    error: str


class AssistantTextBlock(TypedDict):
    type: Literal["text"]
    text: str


class AssistantThinkingBlock(TypedDict):
    type: Literal["thinking"]
    thinking: str


AssistantStreamBlock = AssistantTextBlock | AssistantThinkingBlock


class StreamCompleteErrorPayload(TypedDict):
    message: str
    should_contact_support: bool
    should_clear_history: bool


class ErrorHistoryPayload(TypedDict, total=False):
    message: str
    should_contact_support: bool
    should_clear_history: bool
    raw_error: str
    assistant_error_type: str


AnthropicStreamEvent = (
    RawMessageStartEvent
    | RawMessageDeltaEvent
    | RawMessageStopEvent
    | RawContentBlockStartEvent
    | RawContentBlockDeltaEvent
    | RawContentBlockStopEvent
)


@dataclass
class AgentHarness:
    conn: SocketIoThread
    claude: ClaudeSDKClient | None = None
    system_prompt: str | None = None
    pending_operations: dict[str, asyncio.Future] = field(default_factory=dict)
    executing_cells: set[str] = field(default_factory=set)
    operation_counter: int = 0
    current_request_id: str | None = None
    should_auto_continue: bool = False
    pending_auto_continue: bool = False
    pause_until_user_query: bool = False

    pending_messages: asyncio.Queue = field(default_factory=asyncio.Queue)
    conversation_task: asyncio.Task | None = None
    conversation_running: bool = False
    agent_session_id: int | None = None
    claude_session_id: str | None = None
    latest_notebook_context: dict = field(default_factory=dict)
    current_status: str | None = None
    expected_widgets: dict[str, Any | None] = field(default_factory=dict)
    behavior: Behavior = Behavior.step_by_step
    latest_notebook_state: str | None = None
    current_plan: dict | None = None
    in_memory_history: list[dict] = field(default_factory=list)
    mcp_server: McpSdkServerConfig = field(default_factory=lambda: agent_tools_mcp)
    mcp_allowed_tools: list[str] = field(
        default_factory=lambda: list(MCP_ALLOWED_TOOL_NAMES)
    )
    current_query_task: asyncio.Task | None = None
    open_stream_blocks: dict[int, str] = field(default_factory=dict)

    async def send(self, msg: dict[str, Any]) -> None:
        msg_type = msg.get("type", "unknown")
        if msg_type != "agent_stream_delta":
            print(f"[agent] Sending message: {msg_type}")

        await self.conn.send(msg)

    async def _notify_history_updated(self, *, request_id: str | None = None) -> None:
        await self.send({
            "type": "agent_history_updated",
            "session_id": str(self.agent_session_id)
            if self.agent_session_id is not None
            else None,
            **({"request_id": request_id} if request_id else {}),
        })

    async def _fetch_history_from_db(self) -> list[dict]:
        if skip_db_history:
            return self.in_memory_history
        assert self.agent_session_id is not None
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
        return [
            {
                "payload": n.get("payload"),
                "request_id": n.get("requestId"),
                "template_version_id": n.get("templateVersionId"),
            }
            for n in nodes
        ]

    async def _build_messages_from_db(
        self, *, exclude_request_id: str | None = None
    ) -> list[MessageParam]:
        history = await self._fetch_history_from_db()
        anthropic_messages: list[MessageParam] = []

        for item in history:
            if (
                exclude_request_id is not None
                and item.get("request_id") == exclude_request_id
            ):
                continue

            payload = item.get("payload")

            t = payload.get("type") if isinstance(payload, dict) else None
            if t == "anthropic_message":
                role = payload.get("role")
                content = payload.get("content")

                display_widgets = payload.get("display_widgets")
                if (
                    role == "user"
                    and display_widgets is not None
                    and isinstance(content, str)
                ):
                    for widget in display_widgets:
                        ref_pattern = f"@({widget['widgetKey']}|{widget['id']})"
                        inline_widget = f'<Widget label="{widget["label"]}" type="{widget["widgetType"]}" widget_key="{widget["widgetKey"]}" cell="{widget["cellDisplayName"]}"/>'
                        content = content.replace(ref_pattern, inline_widget)

                if (
                    role == "user"
                    and isinstance(content, dict)
                    and content.get("type") == "cell_result"
                ):
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
                        if (
                            isinstance(block, dict)
                            and block.get("type") == "thinking_summary"
                        ):
                            continue

                        if (
                            isinstance(block, dict)
                            and block.get("type") == "tool_result"
                        ):
                            block = block.copy()
                            block_content = block.get("content", "{}")
                            if isinstance(block_content, str):
                                try:
                                    result = json.loads(block_content)
                                except json.JSONDecodeError:
                                    print(
                                        f"[agent] tool_result parse failed tool_use_id={block.get('tool_use_id', '?')}"
                                    )
                                    result = None
                                if (
                                    isinstance(result, dict)
                                    and "original_code" in result
                                ):
                                    result.pop("original_code")
                                    block["content"] = json.dumps(
                                        result, sort_keys=True
                                    )
                        cleaned_content.append(block)
                    content = cleaned_content

                template_version_id = item.get("template_version_id")
                if role == "user" and template_version_id is not None:
                    checkpoint_content = f"[auto-generated metadata] template_version_id={template_version_id}"
                    anthropic_messages.append({
                        "role": "user",
                        "content": checkpoint_content,
                    })

                if role in {"user", "assistant"} and (isinstance(content, (str, list))):
                    anthropic_messages.append({"role": role, "content": content})

            elif t == "cancellation":
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
                is_tool_result_msg = isinstance(content, list) and any(
                    isinstance(b, dict) and b.get("type") == "tool_result"
                    for b in content
                )
                if is_tool_result_msg:
                    reordered.append(msg)
                    for block in content:
                        if (
                            isinstance(block, dict)
                            and block.get("type") == "tool_result"
                        ):
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

    @staticmethod
    def _truncate_legacy_history_message(msg: MessageParam) -> MessageParam:
        content = msg.get("content")

        if isinstance(content, str):
            if len(content) > 500:
                return {
                    "role": msg["role"],
                    "content": content[:500] + "...[truncated]",
                }
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
                    block_content = block.get("content", "{}")

                    if isinstance(block_content, list):
                        text_blocks = [
                            b
                            for b in block_content
                            if isinstance(b, dict) and b.get("type") == "text"
                        ]
                        if len(text_blocks) > 0:
                            result_str = text_blocks[0].get("text", "{}")
                            result = json.loads(result_str)
                            result["_note"] = "image removed during truncation"
                            truncated_blocks.append({
                                "type": "tool_result",
                                "tool_use_id": block.get("tool_use_id"),
                                "content": json.dumps(result),
                            })
                        else:
                            truncated_blocks.append(block)
                        continue

                    result_str = block_content
                    result = json.loads(result_str)
                    tool_name = result.get("tool_name")

                    if tool_name in {"create_cell", "create_markdown_cell"}:
                        result = {k: v for k, v in result.items() if k != "code"}
                    elif tool_name == "edit_cell":
                        result = {
                            k: v
                            for k, v in result.items()
                            if k not in {"code", "original_code"}
                        }
                    elif tool_name == "bash":
                        result = {k: v for k, v in result.items() if k != "output"}
                    elif tool_name == "execute_code":
                        result = {k: v for k, v in result.items() if k != "code"}
                        if (
                            "stdout" in result
                            and isinstance(result.get("stdout"), str)
                            and len(result["stdout"]) > 1500
                        ):
                            result["stdout"] = (
                                "...[truncated]\n" + result["stdout"][-1500:]
                            )
                        if (
                            "stderr" in result
                            and isinstance(result.get("stderr"), str)
                            and len(result["stderr"]) > 1500
                        ):
                            result["stderr"] = (
                                "...[truncated]\n" + result["stderr"][-1500:]
                            )
                        if (
                            "exception" in result
                            and isinstance(result.get("exception"), str)
                            and len(result["exception"]) > 1500
                        ):
                            result["exception"] = (
                                "...[truncated]\n" + result["exception"][-1500:]
                            )

                    truncated_blocks.append({
                        "type": "tool_result",
                        "tool_use_id": block.get("tool_use_id"),
                        "content": json.dumps(result),
                    })
                    continue

                if block_type == "tool_use":
                    tool_name = block.get("name")
                    inp = block.get("input")
                    truncated_inp = inp

                    if isinstance(inp, dict) and tool_name is not None:
                        if tool_name in {"create_cell", "create_markdown_cell"}:
                            truncated_inp = {
                                k: v for k, v in inp.items() if k != "code"
                            }
                        elif tool_name == "edit_cell":
                            truncated_inp = {
                                k: v for k, v in inp.items() if k != "new_code"
                            }
                        elif tool_name == "update_plan":
                            truncated_inp = {
                                k: v
                                for k, v in inp.items()
                                if k not in {"plan", "plan_diff"}
                            }

                    truncated_blocks.append({
                        "type": "tool_use",
                        "id": block.get("id"),
                        "name": tool_name,
                        "input": truncated_inp,
                    })
                    continue

                if block_type == "text":
                    text = block.get("text", "")
                    if len(text) > 1000:
                        truncated_blocks.append({
                            "type": "text",
                            "text": text[:1000] + "...[truncated]",
                        })
                    else:
                        truncated_blocks.append(block)
                    continue

                truncated_blocks.append(block)

            return {"role": msg["role"], "content": truncated_blocks}

        return msg

    def _build_legacy_history_message(self, msg: MessageParam) -> str:
        role = msg["role"]
        content = msg["content"]

        if isinstance(content, str):
            return f"{role}:\n{content}"

        block_lines: list[str] = []
        for block in content:
            block_type = block.get("type")
            if block_type == "text":
                block_lines.append(block["text"])
                continue

            if block_type == "tool_use":
                tool_name = block.get("name", "unknown")
                tool_input = json.dumps(block.get("input"), sort_keys=True, default=str)
                block_lines.append(f"[tool_use] {tool_name} input={tool_input}")
                continue

            if block_type == "tool_result":
                block_content = block.get("content", "")
                if isinstance(block_content, list):
                    tool_result_text = self._extract_tool_result_text(block_content)
                    content_text = (
                        tool_result_text
                        if tool_result_text is not None
                        else json.dumps(block_content, default=str)
                    )
                else:
                    content_text = str(block_content)
                prefix = (
                    "[tool_result error]" if block.get("is_error") else "[tool_result]"
                )
                block_lines.append(f"{prefix} {content_text}")
                continue

            block_lines.append(json.dumps(block, default=str))

        return f"{role}:\n" + "\n".join(block_lines)

    async def _build_legacy_history_prompt_block(
        self, *, exclude_request_id: str | None
    ) -> str | None:
        if self.claude_session_id is not None:
            return None

        messages = await self._build_messages_from_db(
            exclude_request_id=exclude_request_id
        )
        if len(messages) == 0:
            return None

        windowed_messages = messages[-legacy_history_window_size:]
        truncated_messages = [
            self._truncate_legacy_history_message(msg) for msg in windowed_messages
        ]
        included_messages = [
            self._build_legacy_history_message(msg) for msg in truncated_messages
        ]
        history_text = "\n".join(included_messages)

        print(
            "[agent] Injecting legacy history into first turn "
            f"windowed_messages len={len(windowed_messages)})"
        )

        return (
            "This is prior conversation context. Older turns may be omitted.\n"
            f"{history_text}"
        )

    async def _insert_history(
        self,
        *,
        event_type: str = "anthropic_message",
        role: HistoryRole = "user",
        payload: dict[str, Any],
        request_id: str | None = None,
        tx_id: str | None = None,
        template_version_id: str | None = None,
    ) -> None:
        payload = {
            "type": event_type,
            "role": role,
            "timestamp": int(time.time() * 1000),
            **payload,
        }

        if skip_db_history:
            self.in_memory_history.append({
                "payload": payload,
                "request_id": request_id,
                "template_version_id": template_version_id,
            })
            return

        assert self.agent_session_id is not None

        variables = {
            "sessionId": str(self.agent_session_id),
            "eventType": event_type,
            "payload": payload,
            "requestId": request_id,
            "txId": tx_id,
            "templateVersionId": str(template_version_id)
            if template_version_id is not None
            else None,
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
        self, tool_response: str | list[dict[str, Any]] | None
    ) -> ToolResultContent:
        # todo(rteqs): this is kinda sus. why is there a early return in the for loop.
        if isinstance(tool_response, str):
            return tool_response

        if isinstance(tool_response, list):
            normalized_blocks: list[TextToolResultBlock | ImageToolResultBlock] = []
            for item in tool_response:
                if isinstance(item, dict):
                    block_type = item.get("type")
                    text_value = item.get("text")
                    source = item.get("source")
                    image_data = item.get("data")
                    mime_type = item.get("mimeType")
                else:
                    block_type = getattr(item, "type", None)
                    text_value = getattr(item, "text", None)
                    source = getattr(item, "source", None)
                    image_data = getattr(item, "data", None)
                    mime_type = getattr(item, "mimeType", None)

                if block_type == "text" and isinstance(text_value, str):
                    normalized_blocks.append({"type": "text", "text": text_value})
                    continue

                if block_type == "image":
                    if isinstance(source, dict):
                        source_type = source.get("type")
                        media_type = source.get("media_type")
                        data = source.get("data")
                    else:
                        source_type = getattr(source, "type", None)
                        media_type = getattr(source, "media_type", None)
                        data = getattr(source, "data", None)

                    if (
                        source_type == "base64"
                        and isinstance(media_type, str)
                        and isinstance(data, str)
                    ):
                        normalized_blocks.append({
                            "type": "image",
                            "source": {
                                "type": "base64",
                                "media_type": media_type,
                                "data": data,
                            },
                        })
                        continue

                    if isinstance(mime_type, str) and isinstance(image_data, str):
                        normalized_blocks.append({
                            "type": "image",
                            "source": {
                                "type": "base64",
                                "media_type": mime_type,
                                "data": image_data,
                            },
                        })
                        continue

                return json.dumps(tool_response, default=str)

            if len(normalized_blocks) > 0:
                return normalized_blocks

            return json.dumps(tool_response, default=str)

        return json.dumps(tool_response, default=str)

    def _extract_tool_result_text(self, content: ToolResultContent) -> str | None:
        if isinstance(content, str):
            return content
        for block in content:
            if block["type"] == "text":
                return block["text"]
        return None

    def _build_generic_tool_result_content(
        self, *, tool_name: str, content: ToolResultContent, is_error: bool
    ) -> ToolResultContent:
        text_content = self._extract_tool_result_text(content)
        if text_content is None:
            return content

        try:
            parsed = json.loads(text_content)
        except Exception:
            parsed = None
        if isinstance(parsed, dict) and isinstance(parsed.get("tool_name"), str):
            return content

        result_payload: GenericToolResultPayload = {
            "tool_name": tool_name,
            "success": not is_error,
            "summary": (
                f"{tool_name} failed" if is_error else f"{tool_name} completed"
            ),
        }
        if text_content.strip() != "":
            result_payload["raw_output"] = text_content
        if is_error:
            result_payload["error"] = text_content
            if text_content.strip() != "":
                result_payload["summary"] = f"{tool_name} failed: {text_content[:200]}"

        return [{"type": "text", "text": json.dumps(result_payload)}]

    async def _post_tool_use_attach_notebook_state(
        self,
        input_data: PostToolUseHookInput,
        tool_use_id: str | None,
        _context: HookContext,
    ) -> SyncHookJSONOutput:
        tool_name = input_data["tool_name"]
        notebook_state = await self.refresh_cells_context()
        self.latest_notebook_state = notebook_state

        print(
            "[agent] PostToolUse attached notebook state context "
            f"for {tool_name} (tool_use_id={tool_use_id})"
        )
        return {
            "hookSpecificOutput": {
                "hookEventName": "PostToolUse",
                "additionalContext": (
                    "<current_notebook_state>\n"
                    f"{notebook_state}\n"
                    "</current_notebook_state>"
                ),
            }
        }

    def _extract_tool_blocks_from_sdk_message(
        self, msg: AssistantMessage | UserMessage
    ) -> tuple[list[ToolUseHistoryBlock], list[ToolResultHistoryBlock]]:
        content = msg.content
        if isinstance(content, str):
            return [], []

        tool_use_blocks: list[ToolUseHistoryBlock] = []
        tool_result_blocks: list[ToolResultHistoryBlock] = []

        for block in content:
            if isinstance(block, ToolUseBlock):
                tool_use_id = block.id
                tool_name = block.name
                if tool_use_id == "" or tool_name == "":
                    continue
                tool_use_blocks.append({
                    "type": "tool_use",
                    "id": tool_use_id,
                    "name": tool_name,
                    "input": block.input,
                })
                continue

            if not isinstance(block, ToolResultBlock):
                continue

            tool_use_id = block.tool_use_id
            if tool_use_id == "":
                continue

            normalized_content = self._normalize_tool_result_content(block.content)

            payload_block = {
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
        msg: AssistantMessage | UserMessage,
        request_id: str | None,
        tool_use_index: dict[str, str],
    ) -> None:
        tool_use_blocks, tool_result_blocks = (
            self._extract_tool_blocks_from_sdk_message(msg)
        )
        message_role: Literal["assistant", "user"] = (
            "assistant" if isinstance(msg, AssistantMessage) else "user"
        )

        for block in tool_use_blocks:
            tool_use_index[block["id"]] = block["name"]

        if len(tool_use_blocks) > 0:
            try:
                await self._insert_history(
                    role=message_role,
                    payload={"content": tool_use_blocks},
                    request_id=request_id,
                )
            except Exception as e:
                print(f"[agent] Failed to persist message tool_use blocks: {e!s}")

        if len(tool_result_blocks) > 0:
            normalized_tool_result_blocks: list[ToolResultHistoryBlock] = []
            for block in tool_result_blocks:
                normalized_tool_result_block: ToolResultHistoryBlock = {
                    "type": "tool_result",
                    "tool_use_id": block["tool_use_id"],
                    "content": block["content"],
                }
                if "is_error" in block:
                    normalized_tool_result_block["is_error"] = block["is_error"]

                tool_name = tool_use_index.get(block["tool_use_id"])
                if tool_name is not None:
                    normalized_tool_result_block["content"] = (
                        self._build_generic_tool_result_content(
                            tool_name=tool_name,
                            content=block["content"],
                            is_error=block.get("is_error", False),
                        )
                    )
                normalized_tool_result_blocks.append(normalized_tool_result_block)

            try:
                await self._insert_history(
                    role=message_role,
                    payload={"content": normalized_tool_result_blocks},
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

    def _load_behavior_context(self) -> tuple[str, str]:
        behavior_file = (
            "proactive.md" if self.behavior == Behavior.proactive else "step_by_step.md"
        )
        turn_behavior_content = (
            context_root / "turn_behavior" / behavior_file
        ).read_text()
        examples_content = (context_root / "examples" / behavior_file).read_text()
        return turn_behavior_content, examples_content

    def _compose_turn_system_prompt(self) -> str:
        self.system_prompt = (context_root.parent / "system_prompt.md").read_text()
        assert self.system_prompt is not None

        turn_behavior_content, examples_content = self._load_behavior_context()

        final_system_prompt = self.system_prompt.replace(
            "TURN_BEHAVIOR_PLACEHOLDER",
            f"<turn_behavior>\n{turn_behavior_content}\n</turn_behavior>",
        )
        return final_system_prompt.replace(
            "EXAMPLES_PLACEHOLDER", f"<examples>\n{examples_content}\n</examples>"
        )

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
                "error": f"OPERATION FAILED: '{action}' timed out after {timeout} seconds. This operation did NOT complete.",
                "tx_id": tx_id,
            }
        finally:
            ret = self.pending_operations.pop(tx_id, None)

            if ret is not None:
                duration = time.time() - start_time
                print(f"[agent] {action} took {duration:.3f}s")

    async def handle_action_response(self, msg: dict[str, Any]) -> None:
        tx_id = msg.get("tx_id")
        assert tx_id is not None

        fut = self.pending_operations.get(tx_id)
        if fut and not fut.done():
            fut.set_result(msg)

    async def _send_usage_update(self, usage: Usage | MessageDeltaUsage) -> None:
        await self.send({
            "type": "agent_usage_update",
            "input_tokens": usage.input_tokens if usage.input_tokens is not None else 0,
            "cache_read_input_tokens": usage.cache_read_input_tokens,
            "cache_creation_input_tokens": usage.cache_creation_input_tokens,
            "context_limit": 200_000,
        })

    def _normalize_claude_session_id(self, raw_session_id: str | None) -> str | None:
        if raw_session_id is None or raw_session_id == "":
            return None
        return raw_session_id

    async def _close_open_stream_blocks(self) -> None:
        if len(self.open_stream_blocks) == 0:
            return

        pending_blocks = sorted(self.open_stream_blocks)
        print(
            f"[agent] Closing {len(pending_blocks)} dangling stream block(s): {pending_blocks}"
        )
        for block_index in pending_blocks:
            await self.send({
                "type": "agent_stream_block_stop",
                "block_index": block_index,
            })
        self.open_stream_blocks.clear()

    async def _wait_for_running_query_to_stop(
        self, timeout_seconds: float = 2.0
    ) -> bool:
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
        self.conversation_task = asyncio.create_task(self.run_agent_loop())

        def _task_done_callback(task: asyncio.Task) -> None:
            try:
                task.result()
            except asyncio.CancelledError:
                pass
            except Exception as e:
                print(f"[agent] conversation_task raised exception: {e}")
                traceback.print_exc()

        self.conversation_task.add_done_callback(_task_done_callback)

    async def _wait_for_message(self) -> tuple[str, str | None]:
        print("[agent] _wait_for_message: waiting for message...")
        while True:
            msg = await self.pending_messages.get()
            msg_type = msg.get("type")

            if self.pause_until_user_query and msg_type in {
                "resume",
                "cell_result",
                "set_widget_value",
            }:
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
                self.current_request_id = (
                    request_id if isinstance(request_id, str) else None
                )
                self.current_status = "thinking"
                self.pause_until_user_query = False

                payload = {"content": msg["content"]}

                if msg.get("display_query") is not None:
                    payload["display_query"] = msg["display_query"]
                if msg.get("display_nodes") is not None:
                    payload["display_nodes"] = msg["display_nodes"]
                if msg.get("display_widgets") is not None:
                    payload["display_widgets"] = msg["display_widgets"]
                if msg.get("hidden") is not None:
                    payload["hidden"] = msg["hidden"]

                template_version_id = msg.get("template_version_id")

                await self._insert_history(
                    payload=payload,
                    request_id=self.current_request_id,
                    template_version_id=template_version_id,
                )

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
                    result_message = (
                        f"✓ Cell {cell_name} ({cell_id}) executed successfully"
                    )
                    result_content = {
                        "type": "cell_result",
                        "message": result_message,
                        "cell_id": cell_id,
                        "cell_name": cell_name,
                        "success": success,
                        "logs": logs,
                    }
                    print(f"[agent] Cell {cell_id} succeeded")
                else:
                    exception = msg.get("exception", '"Unknown error"')
                    result_message = f"✗ Cell {cell_name} ({cell_id}) execution failed"
                    result_content = {
                        "type": "cell_result",
                        "message": result_message,
                        "cell_id": cell_id,
                        "cell_name": cell_name,
                        "success": False,
                        "exception": exception,
                        "logs": logs,
                    }
                    print(f"[agent] Cell {cell_id} failed")

                await self._insert_history(payload={"content": result_content})

                prompt_content = result_message
                if not success:
                    exception = result_content.get("exception")
                    if exception:
                        prompt_content = f"{prompt_content}\n\nException: {exception}"
                if logs:
                    prompt_content = f"{prompt_content}\n\nLogs:\n{logs}"

                if self.pending_auto_continue and len(self.executing_cells) == 0:
                    self.pending_auto_continue = False
                    return (
                        "Continue with the next step.\n\n" + prompt_content,
                        self.current_request_id,
                    )

                return (
                    "Cell execution update:\n\n" + prompt_content,
                    self.current_request_id,
                )

            if msg_type == "set_widget_value":
                data = msg.get("data", {})
                widget_info = ", ".join(f"{k}={v}" for k, v in data.items())
                content = f"User provided input via widget(s): {widget_info}"

                await self._insert_history(payload={"content": content, "hidden": True})

                return content, self.current_request_id

            print(f"[agent] Unknown pending message type={msg_type}, ignoring")

    async def _complete_turn(self) -> None:
        if self.current_request_id is None:
            return

        should_continue = self.should_auto_continue
        self.should_auto_continue = False
        self.pending_auto_continue = False

        await self._notify_history_updated(request_id=self.current_request_id)

        if should_continue:
            print("[agent] Auto-continuing as requested by model")
            await self.pending_messages.put({
                "type": "user_query",
                "content": "Continue with the next step.",
                "request_id": self.current_request_id,
                "hidden": True,
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
        if self.claude is None:
            return
        try:
            await self.claude.disconnect()
        except Exception:
            traceback.print_exc()
        finally:
            self.claude = None

    async def _load_claude_session_id(self, session_id: int) -> str | None:
        try:
            response = await gql_query(
                auth=auth_token_sdk,
                query="""
                    query AgentSessionClaudeSessionId($id: BigInt!) {
                        agentSession(id: $id) {
                            claudeSessionId
                        }
                    }
                """,
                variables={"id": session_id},
            )
            parsed_response = validate(response, AgentSessionQueryResp)
            return self._normalize_claude_session_id(
                parsed_response["data"]["agentSession"]["claudeSessionId"]
            )
        except Exception as e:
            print(f"[agent] Failed to load agent claude_session_id: {e!s}")
            return None

    async def _persist_claude_session_id(self, claude_session_id: str) -> None:
        if self.agent_session_id is None:
            return

        await gql_query(
            auth=auth_token_sdk,
            query="""
                mutation UpdateAgentSessionClaudeSessionId($id: BigInt!, $claudeSessionId: UUID) {
                    updateAgentSession(input: {id: $id, patch: {claudeSessionId: $claudeSessionId}}) {
                        clientMutationId
                    }
                }
            """,
            variables={
                "id": self.agent_session_id,
                "claudeSessionId": claude_session_id,
            },
        )

    async def _capture_claude_session_id(self, raw_session_id: str | None) -> None:
        session_id = self._normalize_claude_session_id(raw_session_id)
        if session_id is None or self.agent_session_id is None:
            return

        if self.claude_session_id == session_id:
            return

        self.claude_session_id = session_id
        try:
            await self._persist_claude_session_id(session_id)
        except Exception as e:
            print(f"[agent] Failed to persist Claude session id: {e!s}")
            return

        print(
            "[agent] Persisted Claude session id "
            f"(db_session_id={self.agent_session_id}, claude_session_id={session_id})"
        )

    def _parse_stream_event(
        self, raw_event: dict[str, Any]
    ) -> AnthropicStreamEvent | None:
        event_type = raw_event.get("type")
        try:
            if event_type == "message_start":
                return RawMessageStartEvent(**raw_event)
            if event_type == "message_delta":
                return RawMessageDeltaEvent(**raw_event)
            if event_type == "message_stop":
                return RawMessageStopEvent(**raw_event)
            if event_type == "content_block_start":
                return RawContentBlockStartEvent(**raw_event)
            if event_type == "content_block_delta":
                return RawContentBlockDeltaEvent(**raw_event)
            if event_type == "content_block_stop":
                return RawContentBlockStopEvent(**raw_event)
        except Exception as e:
            print(f"[agent] Invalid stream event payload type={event_type}: {e!s}")
            return None
        print(f"[agent] Unknown stream event type={event_type}")
        return None

    def _build_sdk_env(self) -> dict[str, str]:
        sdk_base_url = f"{nucleus_url}/infer/plots-agent/anthropic"
        sdk_env = {
            "ANTHROPIC_BASE_URL": sdk_base_url,
            "ANTHROPIC_API_KEY": sdk_token,
            "ANTHROPIC_CUSTOM_HEADERS": f"Pod-Id: {pod_id!s}",
        }

        print(f"[agent] SDK gateway mode: nucleus-proxy ({sdk_base_url})")
        return sdk_env

    async def _connect_sdk_client(self, *, resume_session_id: str | None) -> None:
        self.system_prompt = self._compose_turn_system_prompt()
        sdk_env = self._build_sdk_env()

        self.claude = ClaudeSDKClient(
            options=ClaudeAgentOptions(
                system_prompt=self.system_prompt,
                include_partial_messages=True,
                mcp_servers={MCP_SERVER_NAME: self.mcp_server},
                allowed_tools=[*self.mcp_allowed_tools, *SDK_BUILTIN_ALLOWED_TOOLS],
                permission_mode="acceptEdits",
                model="claude-opus-4-6",
                thinking={"type": "adaptive"},
                resume=resume_session_id,
                env=sdk_env,
                hooks={
                    "PostToolUse": [
                        HookMatcher(
                            matcher=NOTEBOOK_MUTATION_TOOL_MATCHER,
                            hooks=[self._post_tool_use_attach_notebook_state],
                        )
                    ]
                },
            )
        )
        await self.claude.connect()

        try:
            mcp_status = await self.claude.get_mcp_status()
            print(f"[agent] MCP status: {json.dumps(mcp_status)}")
        except Exception as e:
            print(f"[agent] Failed to fetch MCP status: {e!s}")

    async def _reset_for_new_session(self) -> None:
        await self._stop_conversation_loop()

        if len(self.pending_operations) > 0:
            print(
                f"[agent] Failing {len(self.pending_operations)} pending operation(s) for new session"
            )
            for future in self.pending_operations.values():
                if not future.done():
                    future.set_result({
                        "status": "error",
                        "error": "Operation interrupted due to session change. Retry in the new session.",
                    })
            self.pending_operations.clear()

        self.pause_until_user_query = False
        self.executing_cells.clear()
        self.expected_widgets.clear()
        self.open_stream_blocks.clear()
        self.current_request_id = None
        self.current_status = None
        self.should_auto_continue = False
        self.pending_auto_continue = False
        self.current_plan = None

        while not self.pending_messages.empty():
            try:
                self.pending_messages.get_nowait()
            except asyncio.QueueEmpty:
                break

        if skip_db_history:
            self.in_memory_history.clear()

        for file in (context_root / "agent_scratch").rglob("*"):
            if file.name == ".gitkeep":
                continue
            file.unlink()

        self.claude_session_id = None

    async def _handle_stream_event(
        self,
        event: AnthropicStreamEvent,
        *,
        collected_assistant_blocks: dict[int, AssistantStreamBlock] | None = None,
    ) -> None:
        if event.type == "message_start":
            self.open_stream_blocks.clear()
            return

        if event.type == "content_block_start":
            block_index = event.index
            block = event.content_block

            if block_index in self.open_stream_blocks:
                await self.send({
                    "type": "agent_stream_block_stop",
                    "block_index": block_index,
                })
            self.open_stream_blocks[block_index] = block.type

            if collected_assistant_blocks is not None:
                if block.type == "text":
                    collected_assistant_blocks[block_index] = {
                        "type": "text",
                        "text": "",
                    }
                elif block.type == "thinking":
                    collected_assistant_blocks[block_index] = {
                        "type": "thinking",
                        "thinking": "",
                    }

            payload: dict[str, Any] = {
                "type": "agent_stream_block_start",
                "block_index": block_index,
                "block_type": block.type,
            }
            if block.type == "tool_use":
                payload["block_id"] = block.id
                payload["block_name"] = block.name
                print(
                    "[agent] tool_use block started "
                    f"index={block_index} name={payload['block_name']} id={payload['block_id']}"
                )

            await self.send(payload)
            return

        if event.type == "content_block_delta":
            block_index = event.index
            delta = event.delta

            if delta.type == "text_delta":
                await self.send({
                    "type": "agent_stream_delta",
                    "block_index": block_index,
                    "block_type": "text",
                    "delta": delta.text,
                })

            elif delta.type == "thinking_delta":
                await self.send({
                    "type": "agent_stream_delta",
                    "block_index": block_index,
                    "block_type": "thinking",
                    "delta": delta.thinking,
                })

            elif delta.type == "input_json_delta":
                await self.send({
                    "type": "agent_stream_delta",
                    "block_index": block_index,
                    "block_type": "tool_use",
                    "delta": delta.partial_json,
                })

            if collected_assistant_blocks is not None:
                existing = collected_assistant_blocks.get(block_index)
                if existing is not None:
                    if delta.type == "text_delta" and existing["type"] == "text":
                        existing["text"] += delta.text
                    elif (
                        delta.type == "thinking_delta"
                        and existing["type"] == "thinking"
                    ):
                        existing["thinking"] += delta.thinking
            return

        if event.type == "content_block_stop":
            block_index = event.index
            if block_index in self.open_stream_blocks:
                self.open_stream_blocks.pop(block_index, None)
            await self.send({
                "type": "agent_stream_block_stop",
                "block_index": block_index,
            })
            return

        if event.type == "message_delta":
            return

        if event.type == "message_stop":
            await self._close_open_stream_blocks()
            return

    async def _run_query(self, *, prompt: str, request_id: str | None) -> None:
        assert self.claude is not None
        assert self.agent_session_id is not None

        self.open_stream_blocks.clear()
        run_started_at = time.perf_counter()
        assistant_blocks_by_index: dict[int, AssistantStreamBlock] = {}
        assistant_message_started_at: float | None = None
        persisted_assistant_message_this_turn = False
        final_query_error: str | None = None
        stream_complete_error: StreamCompleteErrorPayload | None = None
        usage_data: Usage | MessageDeltaUsage | None = None
        tool_use_index: dict[str, str] = {}
        assistant_error_type: str | None = None
        stream_message_open = False
        buffered_tool_messages: list[AssistantMessage | UserMessage] = []
        final_submit_response_reached = False
        dropping_post_submit_stream_message = False

        async def persist_current_assistant_message(
            *, duration_seconds: float | None = None
        ) -> bool:
            if len(assistant_blocks_by_index) == 0:
                return False

            assistant_content = [
                assistant_blocks_by_index[idx]
                for idx in sorted(assistant_blocks_by_index)
            ]

            resolved_duration = duration_seconds
            if resolved_duration is None:
                if assistant_message_started_at is not None:
                    message_duration_start = assistant_message_started_at
                else:
                    message_duration_start = run_started_at
                resolved_duration = max(
                    0.0, time.perf_counter() - message_duration_start
                )

            assistant_message_persisted = False
            try:
                await self._insert_history(
                    role="assistant",
                    payload={
                        "content": assistant_content,
                        "duration": resolved_duration,
                    },
                    request_id=request_id,
                )
                assistant_message_persisted = True
            except Exception as e:
                print(f"[agent] Failed to persist assistant history: {e!s}")
            finally:
                assistant_blocks_by_index.clear()

            return assistant_message_persisted

        # todo(tim): consider replacing buffer system with simpler system
        async def flush_buffered_tool_messages() -> None:
            if len(buffered_tool_messages) == 0:
                return

            buffered_batch = list(buffered_tool_messages)
            buffered_tool_messages.clear()
            for buffered_msg in buffered_batch:
                await self._persist_tool_blocks_from_sdk_message(
                    msg=buffered_msg,
                    request_id=request_id,
                    tool_use_index=tool_use_index,
                )

        await self.send({
            "type": "agent_stream_start",
            "timestamp": int(time.time() * 1000),
        })

        try:
            query_session_id = self.claude_session_id
            if query_session_id is None:
                query_session_id = "default"

            print(
                "[agent] starting SDK query "
                f"(request_id={request_id}, db_session_id={self.agent_session_id}, "
                f"claude_session_id={query_session_id})"
            )
            try:
                await asyncio.wait_for(
                    self.claude.query(prompt=prompt, session_id=query_session_id),
                    timeout=10.0,
                )
            except TimeoutError:
                final_query_error = "Timed out submitting prompt to Claude runtime"
                await self.send({
                    "type": "agent_error",
                    "error": final_query_error,
                    "fatal": False,
                })
                return
            print(f"[agent] SDK query submitted (request_id={request_id})")
            receive_iter = aiter(self.claude.receive_response())
            while True:
                try:
                    msg = await asyncio.wait_for(anext(receive_iter), timeout=90.0)
                except StopAsyncIteration:
                    break
                except TimeoutError:
                    await self.send({
                        "type": "agent_error",
                        "error": "Timed out waiting for model response from Claude runtime",
                        "fatal": False,
                    })
                    final_query_error = (
                        "Timed out waiting for model response from Claude runtime"
                    )
                    break

                if isinstance(msg, StreamEvent):
                    await self._capture_claude_session_id(msg.session_id)
                    event = self._parse_stream_event(msg.event)

                    if event is None:
                        continue

                    if final_submit_response_reached:
                        if event.type == "message_start":
                            dropping_post_submit_stream_message = True
                            assistant_blocks_by_index.clear()
                            assistant_message_started_at = None
                            stream_message_open = False
                            self.open_stream_blocks.clear()
                            continue
                        if dropping_post_submit_stream_message:
                            if event.type == "message_stop":
                                dropping_post_submit_stream_message = False
                                self.open_stream_blocks.clear()
                            continue

                    if event.type == "message_start":
                        stream_message_open = True
                        assistant_blocks_by_index.clear()
                        assistant_message_started_at = time.perf_counter()
                        usage_data = event.message.usage
                    elif event.type == "message_delta":
                        usage_data = event.usage

                    await self._handle_stream_event(
                        event, collected_assistant_blocks=assistant_blocks_by_index
                    )
                    if event.type == "message_stop":
                        message_duration_seconds: float | None = None
                        if assistant_message_started_at is not None:
                            message_duration_seconds = max(
                                0.0, time.perf_counter() - assistant_message_started_at
                            )
                        assistant_message_persisted = (
                            await persist_current_assistant_message(
                                duration_seconds=message_duration_seconds
                            )
                        )
                        if assistant_message_persisted:
                            persisted_assistant_message_this_turn = True
                        assistant_message_started_at = None
                        stream_message_open = False
                        await flush_buffered_tool_messages()
                elif isinstance(msg, ResultMessage):
                    await self._capture_claude_session_id(msg.session_id)
                    print(
                        "[agent] SDK result message "
                        f"subtype={msg.subtype} is_error={msg.is_error} result={msg.result!r}"
                    )
                    if msg.is_error:
                        final_query_error = (
                            msg.result
                            if msg.result is not None
                            else "Claude query failed"
                        )
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
                            agent_error_message = final_query_error
                        await self.send({
                            "type": "agent_error",
                            "error": agent_error_message,
                            "fatal": False,
                        })
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
                        session_id = msg.data.get("session_id")
                        await self._capture_claude_session_id(
                            session_id if isinstance(session_id, str) else None
                        )
                    continue
                elif isinstance(msg, (AssistantMessage, UserMessage)):
                    if isinstance(msg, AssistantMessage) and msg.error is not None:
                        assistant_error_type = msg.error
                        print(f"[agent] assistant message error={msg.error}")
                    if (
                        not final_submit_response_reached
                        and self.pause_until_user_query
                        and self.current_status in {"done", "awaiting_user_response"}
                    ):
                        final_submit_response_reached = True
                        print(
                            "[agent] Final submit_response detected, interrupting SDK query and suppressing follow-up streamed assistant output"
                        )
                        try:
                            await self.claude.interrupt()
                        except Exception as e:
                            print(
                                f"[agent] Failed to interrupt SDK query after final submit_response: {e!s}"
                            )
                    if stream_message_open or len(self.open_stream_blocks) > 0:
                        buffered_tool_messages.append(msg)
                        continue
                    await self._persist_tool_blocks_from_sdk_message(
                        msg=msg, request_id=request_id, tool_use_index=tool_use_index
                    )
                    continue
                else:
                    raise TypeError(f"Unexpected SDK message type: {type(msg)}")

            print(f"[agent] finished SDK query (request_id={request_id})")
        except Exception as e:
            final_query_error = f"Agent SDK query failed: {e!s}"
            await self.send({
                "type": "agent_error",
                "error": final_query_error,
                "fatal": False,
            })
        finally:
            if assistant_message_started_at is not None:
                pending_message_duration_start = assistant_message_started_at
            else:
                pending_message_duration_start = run_started_at
            assistant_message_persisted = await persist_current_assistant_message(
                duration_seconds=max(
                    0.0, time.perf_counter() - pending_message_duration_start
                )
            )
            if assistant_message_persisted:
                persisted_assistant_message_this_turn = True
            await flush_buffered_tool_messages()

            if final_query_error is not None:
                try:
                    if stream_complete_error is not None:
                        error_history_payload: ErrorHistoryPayload = {
                            "message": stream_complete_error["message"],
                            "should_contact_support": stream_complete_error[
                                "should_contact_support"
                            ],
                            "should_clear_history": stream_complete_error[
                                "should_clear_history"
                            ],
                            "raw_error": final_query_error,
                        }
                        if assistant_error_type is not None:
                            error_history_payload["assistant_error_type"] = (
                                assistant_error_type
                            )
                    else:
                        error_history_payload = {
                            "message": final_query_error,
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
            stream_complete_payload: dict[str, Any] = {"type": "agent_stream_complete"}
            if stream_complete_error is not None:
                stream_complete_payload["error"] = stream_complete_error
            await self.send(stream_complete_payload)

    async def _run_query_with_turn_prompt(
        self, *, query: str, request_id: str | None
    ) -> None:
        assert self.claude is not None

        try:
            turn_behavior_content, examples_content = self._load_behavior_context()
            notebook_state = await self.refresh_cells_context()
            self.latest_notebook_state = notebook_state
            legacy_history_prompt_block = await self._build_legacy_history_prompt_block(
                exclude_request_id=request_id
            )

            context_blocks = [
                f"<turn_behavior>\n{turn_behavior_content}\n</turn_behavior>",
                f"<examples>\n{examples_content}\n</examples>",
            ]
            if legacy_history_prompt_block is not None:
                context_blocks.append(
                    "<legacy_conversation_history>\n"
                    f"{legacy_history_prompt_block}\n"
                    "</legacy_conversation_history>"
                )
            context_blocks.append(
                f"<current_notebook_state>\n{notebook_state}\n</current_notebook_state>"
            )
            if self.current_plan is not None:
                plan_content = json.dumps(self.current_plan, indent=2)
                context_blocks.append(
                    f"<current_plan>\n{plan_content}\n</current_plan>"
                )
            context_blocks.append(f"<user_request>\n{query}\n</user_request>")
            turn_prompt = "\n\n".join(context_blocks)

            print(
                "[agent] Turn prompt ready "
                f"(behavior={self.behavior.value}, "
                f"has_plan={self.current_plan is not None}, "
                f"legacy_history={legacy_history_prompt_block is not None}, "
                f"notebook_chars={len(self.latest_notebook_state or '')})"
            )
        except Exception as e:
            await self.send({
                "type": "agent_error",
                "error": f"Failed to build turn prompt: {e!s}",
                "fatal": False,
            })
            return

        await self._run_query(prompt=turn_prompt, request_id=request_id)

    async def handle_init(self, msg: dict[str, Any]) -> None:
        try:
            new_session_id = int(msg.get("session_id"))
        except Exception:
            print("[agent] invalid session id")
            await self.send({
                "type": "agent_error",
                "error": "Invalid agent session id.",
                "fatal": False,
            })
            return

        raw_init_claude_session_id = msg.get("claude_session_id")
        init_claude_session_id = self._normalize_claude_session_id(
            raw_init_claude_session_id
            if isinstance(raw_init_claude_session_id, str)
            else None
        )

        session_changed = new_session_id != self.agent_session_id
        if session_changed:
            print(f"[agent] Session initialized/changed: {new_session_id}")
            self.agent_session_id = new_session_id

        if self.claude is not None and not session_changed:
            print("[agent] SDK client already initialized; skipping re-init")
            if len(self.pending_operations) > 0:
                print(
                    f"[agent] Failing {len(self.pending_operations)} pending operations for retry"
                )
                for future in self.pending_operations.values():
                    if not future.done():
                        future.set_result({
                            "status": "error",
                            "error": "Connection was reset during operation. Please retry.",
                        })
                self.pending_operations.clear()
            await self.send({"type": "agent_status", "status": "ready"})
            if self.conversation_task is None or self.conversation_task.done():
                self._start_conversation_loop()
            return

        if self.claude is not None and session_changed:
            if (
                self.current_query_task is not None
                and not self.current_query_task.done()
            ):
                print("[agent] Interrupting running query due to session change")
                await self.claude.interrupt()
                await self._wait_for_running_query_to_stop()
            await self._reset_for_new_session()
            await self._disconnect_sdk_client()

        resume_session_id = init_claude_session_id
        if resume_session_id is None:
            resume_session_id = await self._load_claude_session_id(
                self.agent_session_id
            )
        self.claude_session_id = resume_session_id

        if self.claude is None:
            await self._connect_sdk_client(resume_session_id=resume_session_id)
            if resume_session_id is None:
                print("[agent] SDK initialized without resume session")
            else:
                print(
                    f"[agent] SDK initialized with resume session {resume_session_id}"
                )

        await self.send({"type": "agent_status", "status": "ready"})
        print(
            "[agent] Initialization complete "
            f"(db_session_id={self.agent_session_id}, "
            f"claude_session_id={resume_session_id})"
        )
        self._start_conversation_loop()

    async def handle_query(self, msg: dict[str, Any]) -> None:
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
            full_query = f"{query}\n\nHere is the context of the selected nodes the user would like to use: <ContextualNodeData>{json.dumps(contextual_node_data)}</ContextualNodeData>"
        await self.pending_messages.put({
            "type": "user_query",
            "content": full_query,
            "request_id": request_id,
            "display_query": query,
            "display_nodes": contextual_node_data,
            "display_widgets": selected_widgets,
            "template_version_id": template_version_id,
        })

        print(f"[agent] Message queued successfully (request_id={request_id})")

    async def handle_cancel(self, msg: dict[str, Any]) -> None:
        request_id = msg.get("request_id", "unknown")
        print(f"[agent] Cancelling request {request_id}")

        has_running_query = (
            self.current_query_task is not None and not self.current_query_task.done()
        )
        if self.claude is not None and has_running_query:
            await self.claude.interrupt()
            completed_gracefully = await self._wait_for_running_query_to_stop()
            if completed_gracefully:
                print(f"[agent] Cancel completed for request {request_id}")
            else:
                print(
                    f"[agent] Cancel forced task cancellation for request {request_id}"
                )
        else:
            print(f"[agent] No active query to interrupt for request {request_id}")

        if self.current_query_task is not None and self.current_query_task.done():
            self.current_query_task = None

        self.should_auto_continue = False
        self.pending_auto_continue = False
        self.pause_until_user_query = True
        self.current_request_id = None
        self.expected_widgets.clear()
        self.current_status = None
        self.executing_cells.clear()

        if len(self.pending_operations) > 0:
            print(
                f"[agent] Cancelling {len(self.pending_operations)} pending operations"
            )
            for tx_id, future in list(self.pending_operations.items()):
                if not future.done():
                    future.cancel()
                    print(f"[agent]   Cancelled: {tx_id}")
            self.pending_operations.clear()

        while not self.pending_messages.empty():
            try:
                self.pending_messages.get_nowait()
            except asyncio.QueueEmpty:
                break

        for file in (context_root / "agent_scratch").rglob("*"):
            if file.name == ".gitkeep":
                continue
            file.unlink()

        try:
            await self._insert_history(
                event_type="cancellation",
                role="system",
                payload={"reason": "Request cancelled by user"},
                request_id=str(request_id),
            )
        except Exception as e:
            print(f"[agent] Failed to persist cancellation history: {e!s}")

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

    async def update_system_prompt(self, msg: dict[str, Any]) -> dict:
        assert self.claude is not None
        new_content = msg.get("content")
        if not isinstance(new_content, str):
            return {"status": "error", "error": "Invalid content"}

        system_prompt_path = context_root.parent / "system_prompt.md"
        system_prompt_path.write_text(new_content)
        self.system_prompt = new_content
        self.claude.options.system_prompt = new_content

        full_prompt = await self.get_full_prompt()
        return {"status": "success", **full_prompt}

    async def accept(self, msg) -> None:
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

        async with asyncio.TaskGroup() as tg:
            while True:
                try:
                    msg = await harness.conn.recv()
                    tg.create_task(harness.accept(msg))
                except Exception:
                    traceback.print_exc()
                    continue

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
