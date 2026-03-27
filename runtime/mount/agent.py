import asyncio
import json
import os
import socket
import sys
import time
import traceback
import uuid
from collections.abc import Iterable
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Literal, NotRequired, TypedDict

from anthropic.types import (
    MessageDeltaUsage,
    RawContentBlockDeltaEvent,
    RawContentBlockStartEvent,
    RawContentBlockStopEvent,
    RawMessageDeltaEvent,
    RawMessageStartEvent,
    RawMessageStopEvent,
    Usage,
)
from claude_agent_sdk import (
    ClaudeAgentOptions,
    ClaudeSDKClient,
    HookMatcher,
    PermissionResultAllow,
    ResultMessage,
    SystemMessage,
    TextBlock,
    ThinkingBlock,
    ToolResultBlock,
    ToolUseBlock,
    get_session_messages,
)
from claude_agent_sdk.types import (
    AssistantMessage,
    AssistantMessageError,
    HookContext,
    McpSdkServerConfig,
    PostToolUseHookInput,
    PreCompactHookInput,
    StreamEvent,
    SyncHookJSONOutput,
    SystemPromptPreset,
    UserMessage,
)
from latch_data_validation.data_validation import validate
from lplots import _inject
from socketio_thread import SocketIoThread
from tools import MCP_ALLOWED_TOOL_NAMES, MCP_SERVER_NAME, agent_tools_mcp, can_use_tool
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
    "AskUserQuestion",
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

# todo(rteqs): "plan", "build"
Behavior = Literal["step_by_step", "proactive"]

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
    assistant_error_type: str | None
    subtype: str


AnthropicStreamEvent = (
    RawMessageStartEvent
    | RawMessageDeltaEvent
    | RawMessageStopEvent
    | RawContentBlockStartEvent
    | RawContentBlockDeltaEvent
    | RawContentBlockStopEvent
)


class AgentQuery(TypedDict):
    type: Literal["agent_query"]
    request_id: NotRequired[str]
    query: str
    behavior: Behavior
    template_version_id: str | None
    contextual_node_data: NotRequired[dict[str, Any]]  # todo(rteqs): type properly
    selected_widgets: NotRequired[list[dict[str, Any]]]  # todo(rteqs): type properly


class CacheCreation(TypedDict):
    ephemeral_1h_input_tokens: int
    ephemeral_5m_input_tokens: int


class ServerToolUsage(TypedDict):
    web_search_requests: int
    web_fetch_requests: int


@dataclass(frozen=True)
class ResultMessageUsage:
    cache_creation: CacheCreation | None
    cache_creation_input_tokens: int | None
    cache_read_input_tokens: int | None
    input_tokens: int
    output_tokens: int
    server_tool_use: ServerToolUsage | None
    service_tier: Literal["standard", "priority", "batch"] | None
    inference_geo: str
    iterations: list
    speed: str


@dataclass(frozen=True, kw_only=True)
class PendingOperation:
    fut: asyncio.Future
    start_time: float
    action: str
    params: dict


@dataclass
class AgentHarness:
    conn: SocketIoThread
    claude: ClaudeSDKClient | None = None
    system_prompy: str | None = None
    pending_operations: dict[str, PendingOperation] = field(default_factory=dict)
    pending_cells: set[str] = field(default_factory=set)
    pending_widgets: dict[str, Any | None] = field(default_factory=dict)
    pending_question_tx_id: str | None = None
    pending_question_event: asyncio.Event = field(default_factory=asyncio.Event)
    operation_counter: int = 0
    current_request_id: str | None = None

    conversation_task: asyncio.Task | None = None
    agent_session_id: int | None = None  # db_session_id
    claude_session_id: str | None = None
    _needs_reconnect: bool = False
    latest_notebook_context: dict = field(default_factory=dict)
    current_status: str | None = None
    behavior: Behavior = "step_by_step"
    latest_notebook_state: str | None = None
    _last_sent_behavior: Behavior | None = None
    current_plan: dict | None = None
    in_memory_history: list[dict] = field(default_factory=list)
    mcp_server: McpSdkServerConfig = field(default_factory=lambda: agent_tools_mcp)
    mcp_allowed_tools: list[str] = field(
        default_factory=lambda: list(MCP_ALLOWED_TOOL_NAMES)
    )

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
            f"for {tool_name} ({tool_use_id=})"
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

    async def _pre_compact_hook(
        self,
        input_data: PreCompactHookInput,
        tool_use_id: str | None,
        _context: HookContext,
    ) -> SyncHookJSONOutput:
        trigger = input_data.get("trigger", "unknown")

        print(f"[agent] PreCompact hook fired ({trigger=}) ({tool_use_id=}) ")

        await self._insert_history(
            event_type="compaction_start",
            role="system",
            payload=dict(input_data),
            request_id=self.current_request_id,
        )

        return {"continue_": True}

    # todo(rteqs): refactor _normalize_tool_content, _extract_tool_result_text, _build_generic_tool_result_content
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
                else:
                    block_type = getattr(item, "type", None)
                    text_value = getattr(item, "text", None)
                    source = item.get("source")
                    if block_type is None and isinstance(text_value, str):
                        block_type = "text"

                if block_type == "text" and isinstance(text_value, str):
                    normalized_blocks.append({"type": "text", "text": text_value})
                    continue

                if block_type == "image" and isinstance(source, dict):
                    # note(tim): SDK gives image tool res back nested under `source`.
                    media_type = source.get("media_type")
                    data = source.get("data")
                    if (
                        source.get("type") == "base64"
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

    async def create_atomic_operation(
        self, action: str, params: dict | None = None
    ) -> str:
        if params is None:
            params: dict = {}

        self.operation_counter += 1

        tx_id = f"tx_{uuid.uuid4().hex[:12]}"
        loop = asyncio.get_running_loop()
        response_future = loop.create_future()

        start_time = time.time()
        self.pending_operations[tx_id] = PendingOperation(
            fut=response_future, start_time=start_time, action=action, params=params
        )

        try:
            print(f"[agent] -> {action}")
            await self.send({
                "type": "agent_action",
                "action": action,
                "params": params,
                "tx_id": tx_id,
            })
        except Exception:
            self.pending_operations.pop(tx_id, None)
            raise

        return tx_id

    async def wait_for_operation(
        self, tx_id: str, timeout: float | None = 10.0
    ) -> dict:
        assert tx_id in self.pending_operations

        op = self.pending_operations[tx_id]
        force_backend_browser_retry = op.params.get(
            "force_backend_browser_retry", False
        )

        try:
            return await asyncio.wait_for(op.fut, timeout=timeout)
        except asyncio.CancelledError:
            print(
                f"[agent] Operation cancelled (session reinitialized): action={op.action}, tx_id={tx_id}"
            )
            return {
                "status": "error",
                "error": f"OPERATION FAILED: '{op.action}' was interrupted because the session was reinitialized. This operation did NOT complete. You must retry or inform the user.",
            }
        except TimeoutError:
            self.pending_operations.pop(tx_id, None)
            duration = time.time() - op.start_time
            print(f"[agent] {op.action} timed out after {duration:.3f}s")

            if not force_backend_browser_retry:
                print(f"[agent] Retrying {op.action} with backend browser")
                return await self.atomic_operation(
                    op.action,
                    {**op.params, "force_backend_browser_retry": True},
                    timeout=timeout,
                )

            return {
                "status": "error",
                "error": f"OPERATION FAILED: '{op.action}' timed out after {timeout} seconds. This operation did NOT complete.",
                "tx_id": tx_id,
            }
        finally:
            ret = self.pending_operations.pop(tx_id, None)

            if ret is not None:
                duration = time.time() - op.start_time
                print(f"[agent] {op.action} took {duration:.3f}s")

    async def atomic_operation(
        self, action: str, params: dict | None = None, timeout: float | None = 10.0
    ) -> dict:
        try:
            tx_id = await self.create_atomic_operation(action, params)
        except Exception as e:
            return {"status": "error", "error": f"Send failed: {e!s}"}

        return await self.wait_for_operation(tx_id, timeout)

    async def handle_action_response(self, msg: dict[str, Any]) -> None:
        tx_id = msg.get("tx_id")
        assert tx_id is not None

        op = self.pending_operations.get(tx_id)
        if op is None:
            if msg.get("action") == "ask_user_question":
                self.pending_question_event.clear()
                if msg.get("skip") is True:
                    return

                answers = msg.get("answers", {})
                qa_content = {
                    "type": "answers",
                    "content": [
                        {"question": q, "answer": a} for q, a in answers.items()
                    ],
                }
                await self._insert_history(payload={"content": qa_content})
                await self.query(
                    msg=AgentQuery(
                        type="agent_query",
                        behavior=self.behavior,
                        query=json.dumps(
                            asdict(PermissionResultAllow(updated_input=answers))
                        ),
                        template_version_id=None,
                    )
                )

            return

        fut = op.fut
        if fut and not fut.done():
            fut.set_result(msg)

    async def _send_usage_update(
        self, usage: Usage | MessageDeltaUsage | ResultMessageUsage
    ) -> None:
        await self.send({
            "type": "agent_usage_update",
            "input_tokens": usage.input_tokens if usage.input_tokens is not None else 0,
            "cache_read_input_tokens": usage.cache_read_input_tokens,
            "cache_creation_input_tokens": usage.cache_creation_input_tokens,
            "context_limit": 1_000_000,
        })

    async def set_agent_status(self, status: str) -> None:
        self.current_status = status
        await self.send({"type": "agent_status", "status": status})

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
            data = validate(response, AgentSessionQueryResp)["data"]

            if data is None:
                return None

            agent_session = data["agentSession"]
            if agent_session is None:
                return None

            return agent_session["claudeSessionId"]
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

    def _parse_stream_event(
        self, raw_event: dict[str, Any]
    ) -> AnthropicStreamEvent | None:
        # todo(rteqs): use latch_validation
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

    async def connect(self, *, resume_session_id: str | None) -> None:
        self.system_prompt = (context_root.parent / "system_prompt.md").read_text()

        nucleus_llm_url = f"{nucleus_url}/infer/plots-agent/anthropic"
        sdk_env = {
            "ANTHROPIC_BASE_URL": nucleus_llm_url,
            # todo(rteqs): sdk_token here is latch_sdk_token. refactor this so it doesn't confuse
            "ANTHROPIC_API_KEY": sdk_token,
            "ANTHROPIC_CUSTOM_HEADERS": f"Pod-Id: {pod_id!s}",
        }

        self.claude = ClaudeSDKClient(
            options=ClaudeAgentOptions(
                system_prompt=SystemPromptPreset(
                    type="preset", preset="claude_code", append=self.system_prompt
                ),
                include_partial_messages=True,
                mcp_servers={MCP_SERVER_NAME: self.mcp_server},
                allowed_tools=[*self.mcp_allowed_tools, *SDK_BUILTIN_ALLOWED_TOOLS],
                permission_mode="acceptEdits",
                model="claude-opus-4-6",
                thinking={"type": "adaptive"},
                setting_sources=["project"],
                resume=resume_session_id,
                env=sdk_env,
                hooks={
                    "PostToolUse": [
                        HookMatcher(
                            matcher=NOTEBOOK_MUTATION_TOOL_MATCHER,
                            hooks=[self._post_tool_use_attach_notebook_state],
                        )
                    ],
                    "PreCompact": [HookMatcher(hooks=[self._pre_compact_hook])],
                },
                can_use_tool=can_use_tool,
                max_buffer_size=10 * 1024 * 1024,
            )
        )

        # todo(rteqs): validate claude_session_id matches. currently the only way to receive claude_session_id is from messages after making a query.
        # reverse engineer claude code and write to the process directly
        await self.claude.connect()

        try:
            mcp_status = await self.claude.get_mcp_status()
            print(f"[agent] MCP status: {json.dumps(mcp_status)}")
        except Exception as e:
            print(f"[agent] Failed to fetch MCP status: {e!s}")

    async def disconnect(self) -> None:
        if self.claude is None:
            return
        try:
            await self.claude.disconnect()
        except Exception:
            traceback.print_exc()
        finally:
            self.claude = None

    async def handle_init(self, msg: dict[str, Any]) -> None:
        # todo(rteqs): for old sessions that need to be migrated. we need to inject out session history from vacuole into the .jsonl history file
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

        session_changed = new_session_id != self.agent_session_id
        if session_changed:
            print(f"[agent] Session initialized/changed: {new_session_id}")
            self.agent_session_id = new_session_id

        if self.claude is not None and not session_changed:
            print("[agent] SDK client already initialized; skipping re-init")
            await self.send({"type": "agent_status", "status": "ready"})
            return

        if self.claude is not None and session_changed:
            print("[agent] Session changed, will reconnect on next query")
            self._needs_reconnect = True

        resume_session_id = msg.get("claude_session_id")
        if resume_session_id is None or resume_session_id == "":
            resume_session_id = await self._load_claude_session_id(
                self.agent_session_id
            )
        self.claude_session_id = resume_session_id

        await self.send({"type": "agent_status", "status": "ready"})
        print(
            "[agent] Initialization complete "
            f"(db_session_id={self.agent_session_id}, "
            f"claude_session_id={resume_session_id})"
        )

    async def _reset_for_new_session(self) -> None:
        if len(self.pending_operations) > 0:
            print(
                f"[agent] Failing {len(self.pending_operations)} pending operation(s) for new session"
            )
            for op in self.pending_operations.values():
                future = op.fut
                if not future.done():
                    future.set_result({
                        "status": "error",
                        "error": "Operation interrupted due to session change. Retry in the new session.",
                    })
            self.pending_operations.clear()

        self.pending_cells.clear()
        self.pending_widgets.clear()
        self.current_request_id = None
        self.current_status = None
        self.current_plan = None
        self._last_sent_behavior = None

        if skip_db_history:
            self.in_memory_history.clear()

        for file in (context_root / "agent_scratch").rglob("*"):
            if file.name == ".gitkeep":
                continue
            file.unlink()

        self.claude_session_id = None

    async def _handle_stream_event(self, msg: StreamEvent) -> None:
        # note: streaming is just a ui optimization. message delivery here is intentionally unreliable to simplify the backend.
        # client implements logic to stitch streamed blocks and can choose to ignore them altogether.
        # final state of message history is dictated by AssistantMessage and UserMessages

        if self.claude_session_id != msg.session_id:
            return

        event = self._parse_stream_event(msg.event)

        if event is None:
            return

        if event.type == "message_start":
            await self.send({
                "type": "agent_stream_start",
                "timestamp": int(time.time() * 1000),
            })
            return

        if event.type == "message_delta":
            return

        if event.type == "message_stop":
            # todo(rteqs): send message to close stream on the browser
            return

        if event.type == "content_block_start":
            block_index = event.index
            block = event.content_block

            payload = {
                "type": "agent_stream_block_start",
                "block_index": block_index,
                "block_type": block.type,
            }
            if block.type == "tool_use":
                payload |= {"block_id": block.id, "block_name": block.name}
                payload["block_id"] = block.id
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

            return

        if event.type == "content_block_stop":
            await self.send({
                "type": "agent_stream_block_stop",
                "block_index": event.index,
            })
            return

    # todo(rteqs): remove this once we move to "plan", "build"
    def _load_behavior_context(self) -> tuple[str, str]:
        behavior_file = (
            "proactive.md" if self.behavior == "proactive" else "step_by_step.md"
        )
        turn_behavior_content = (
            context_root / "turn_behavior" / behavior_file
        ).read_text()
        examples_content = (context_root / "examples" / behavior_file).read_text()
        return turn_behavior_content, examples_content

    async def create_prompt(self, query: AgentQuery) -> str:
        # todo(rteqs): we should probably move notebook state to a tool or pull it in everytime the model needs it instead of just the prompt
        self.latest_notebook_state = await self.refresh_cells_context()

        context_blocks = []
        if self._last_sent_behavior != self.behavior:
            turn_behavior_content, examples_content = self._load_behavior_context()
            context_blocks.extend([
                f"<turn_behavior>\n{turn_behavior_content}\n</turn_behavior>",
                f"<examples>\n{examples_content}\n</examples>",
            ])
            self._last_sent_behavior = self.behavior

        context_blocks.append(
            f"<current_notebook_state>\n{self.latest_notebook_state}\n</current_notebook_state>"
        )

        if self.current_plan is not None:
            plan_content = json.dumps(self.current_plan, indent=2)
            context_blocks.append(f"<current_plan>\n{plan_content}\n</current_plan>")

        context_blocks.append(f"<user_request>\n{query['query']}\n</user_request>")

        contextual_node_data = query.get("contextual_node_data")
        if contextual_node_data is not None:
            context_blocks.append(
                f"Here is the context of the selected nodes the user would like to use: <contextual_node_data>{json.dumps(contextual_node_data)}</contextual_node_data>"
            )

        selected_widgets = query.get("selected_widgets")
        if selected_widgets is not None:
            context_blocks.append(
                f"Here is the selected widgets the user would like t use: <widgets>{json.dumps(selected_widgets)}<widgets/>"
            )

        return "\n\n".join(context_blocks)

    async def query(self, msg: AgentQuery) -> None:
        if self.claude is not None and self._needs_reconnect:
            print("[agent] Reconnecting SDK for new session")
            try:
                await self.claude.interrupt()
            except Exception:
                pass
            await self._reset_for_new_session()
            await self.disconnect()

        if self.claude is None:
            try:
                await self.connect(resume_session_id=self.claude_session_id)
            except Exception as e:
                print(f"[agent] SDK connection failed: {e!s}")
                traceback.print_exc()
                self.claude = None
                await self.send({
                    "type": "agent_error",
                    "error": f"Agent failed to initialize: {e!s}",
                    "fatal": True,
                })
                return

        self._needs_reconnect = False

        self.behavior = msg["behavior"]
        await self.set_agent_status("thinking")

        prompt = await self.create_prompt(msg)
        request_id = msg.get("request_id")

        print(
            "[agent] starting SDK query "
            f"(request_id={request_id}, db_session_id={self.agent_session_id}, "
            f"claude_session_id={self.claude_session_id})"
        )

        if self.claude_session_id is None:
            await self.claude.query(prompt=prompt)
        else:
            await self.claude.query(prompt=prompt, session_id=self.claude_session_id)

        tool_name_by_tool_use_id: dict[str, str] = {}

        error_type: AssistantMessageError | None = None
        assistant_message_started_at: float | None = None
        try:
            async for res in self.claude.receive_response():
                self.current_request_id = request_id
                if (
                    isinstance(res, SystemMessage)
                    and res.subtype == "init"
                    and self.claude_session_id is None
                    # note(rteqs): this is always the first message for every query
                ):
                    self.claude_session_id = s_id = res.data.get("session_id")

                    if s_id is not None:
                        await self._persist_claude_session_id(s_id)

                if isinstance(res, StreamEvent):
                    event = self._parse_stream_event(res.event)
                    if (
                        event is not None
                        and event.type == "message_start"
                        and assistant_message_started_at is None
                    ):
                        assistant_message_started_at = time.perf_counter()
                    await self._handle_stream_event(res)

                elif isinstance(res, AssistantMessage):
                    turn_duration: float | None = None
                    if assistant_message_started_at is not None:
                        turn_duration = max(
                            0.0, time.perf_counter() - assistant_message_started_at
                        )
                    assistant_message_started_at = None

                    if res.error is not None:
                        error_type = res.error
                        print(f"[agent] assistant message error={res.error}")

                    for c in res.content:
                        if isinstance(c, TextBlock):
                            await self._insert_history(
                                role="assistant",
                                request_id=request_id,
                                payload={
                                    "content": [{"type": "text", "text": c.text}],
                                    **(
                                        {"duration": turn_duration}
                                        if turn_duration is not None
                                        else {}
                                    ),
                                },
                            )

                        elif isinstance(c, ThinkingBlock):
                            await self._insert_history(
                                role="assistant",
                                request_id=request_id,
                                payload={
                                    "content": [
                                        {"type": "thinking", "thinking": c.thinking}
                                    ],
                                    **(
                                        {"duration": turn_duration}
                                        if turn_duration is not None
                                        else {}
                                    ),
                                },
                            )

                        elif isinstance(c, ToolUseBlock):
                            tool_name_by_tool_use_id[c.id] = c.name
                            if c.name == "AskUserQuestion":
                                try:
                                    tx_id = await self.create_atomic_operation(
                                        "ask_user_question",
                                        {"questions": c.input.get("questions", [])},
                                    )
                                    self.pending_question_tx_id = c.input["tx_id"] = (
                                        tx_id
                                    )
                                    self.pending_question_event.set()
                                except Exception:
                                    ...

                            await self._insert_history(
                                role="assistant",
                                request_id=request_id,
                                payload={
                                    "content": [
                                        {
                                            "type": "tool_use",
                                            "id": c.id,
                                            "name": c.name,
                                            "input": c.input,
                                        }
                                    ],
                                    **(
                                        {"duration": turn_duration}
                                        if turn_duration is not None
                                        else {}
                                    ),
                                },
                            )

                elif isinstance(res, UserMessage):
                    if isinstance(res.content, str):
                        await self._insert_history(
                            role="user",
                            request_id=request_id,
                            payload={"content": res.content},
                        )
                        continue

                    for c in res.content:
                        if isinstance(c, TextBlock):
                            await self._insert_history(
                                role="user",
                                request_id=request_id,
                                payload={"content": [{"type": "text", "text": c.text}]},
                            )

                        if isinstance(c, ToolResultBlock):
                            tool_name = tool_name_by_tool_use_id.get(
                                c.tool_use_id, "Unknown"
                            )
                            is_error = c.is_error if c.is_error is not None else False

                            await self._insert_history(
                                role="user",
                                request_id=request_id,
                                payload={
                                    "content": [
                                        {
                                            "type": "tool_result",
                                            "tool_use_id": c.tool_use_id,
                                            "content": self._build_generic_tool_result_content(
                                                tool_name=tool_name,
                                                content=self._normalize_tool_result_content(
                                                    c.content
                                                ),
                                                is_error=is_error,
                                            ),
                                        }
                                    ]
                                },
                            )

                elif isinstance(res, ResultMessage):
                    if res.session_id != self.claude_session_id:
                        continue

                    usage = validate(res.usage, ResultMessageUsage)
                    await self._send_usage_update(usage)

                    if res.subtype == "success":
                        break

                    # todo(rteqs): the old error handling is erroneous. we know from the docs that res.result will be None if subtype != success.
                    # subtype and assistant_error_type (from latest AssistantMessage) alone are too broad to pinpoint what the exact error is and
                    # there doesn't seem to be a fixed schema for res.stop_reason. The below is just a best effort implementation.
                    stop_reason = (
                        res.stop_reason
                        if res.stop_reason is not None
                        else f"unknown error (subtype={res.subtype})"
                    )
                    error_history_payload: ErrorHistoryPayload = {
                        "message": stop_reason,
                        "raw_error": stop_reason,
                        "assistant_error_type": error_type,
                        "subtype": res.subtype,
                        "should_clear_history": "prompt is too long" in stop_reason,
                        "should_contact_support": res.subtype
                        == "error_during_execution",
                    }

                    await self._insert_history(
                        role="system",
                        event_type="error",
                        payload=error_history_payload,
                        request_id=request_id,
                    )
                    break
        except Exception as e:
            print(f"[agent] query failed with unhandled error: {e!s}")
            traceback.print_exc()
            await self.send({
                "type": "agent_error",
                "error": f"Agent query failed: {e!s}",
                "fatal": False,
            })
            try:
                await self._insert_history(
                    event_type="error",
                    role="system",
                    payload={
                        "message": f"Agent query failed: {e!s}",
                        "should_contact_support": True,
                    },
                    request_id=request_id,
                )
            except Exception:
                print("[agent] Failed to persist error history")

    async def interrupt(self, msg: dict[str, Any]) -> None:
        # todo(rteqs): clean up
        request_id = msg.get("request_id", "unknown")
        print(f"[agent] Cancelling request {request_id}")

        if self.claude is not None:
            await self.claude.interrupt()
        else:
            print(f"[agent] No active query to interrupt for request {request_id}")

        self.current_request_id = None
        self.pending_widgets.clear()
        self.current_status = None
        self.pending_cells.clear()

        if len(self.pending_operations) > 0:
            print(
                f"[agent] Cancelling {len(self.pending_operations)} pending operations"
            )
            for tx_id, op in list(self.pending_operations.items()):
                future = op.fut
                if not future.done():
                    future.cancel()
                    print(f"[agent]   Cancelled: {tx_id}")
            self.pending_operations.clear()

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

        messages = []
        if self.claude_session_id is not None:
            messages = get_session_messages(self.claude_session_id)

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
        self.system_prompt = (context_root.parent / "system_prompt.md").read_text()
        self.claude.options.system_prompt = SystemPromptPreset(
            type="preset", preset="claude_code", append=self.system_prompt
        )

        full_prompt = await self.get_full_prompt()
        return {"status": "success", **full_prompt}

    async def accept(self, msg: Any) -> None:
        msg_type = msg.get("type")
        msg_request_id = msg.get("request_id")
        print(
            f"[agent] accept: received message type={msg_type} (request_id={msg_request_id})"
        )

        if msg_type == "init":
            print(f"[agent] Message: {msg_type}")
            await self.handle_init(msg)
        elif msg_type == "agent_query":
            if self.agent_session_id is None:
                print("[agent] Query received before init, ignoring")
                await self.send({
                    "type": "agent_error",
                    "error": "Agent not initialized yet.",
                    "fatal": False,
                })
                return

            query = msg.get("query", "")
            query_preview = query[:60] + "…" if len(query) > 60 else query
            request_id = msg.get("request_id", "unknown")

            print(
                f"[agent] accept: dispatching to handle_query (query={query_preview}, request_id={request_id})"
            )
            await self._insert_history(
                role="user",
                request_id=msg["request_id"],
                payload={"content": msg["query"], "display_query": msg["query"]},
                template_version_id=msg["template_version_id"],
            )
            handle_start = time.time()
            await self.query(msg)
            handle_elapsed = time.time() - handle_start
            print(
                f"[agent] accept: handle_query completed in {handle_elapsed:.3f}s (request_id={request_id})"
            )
        elif msg_type == "agent_cancel":
            request_id = msg.get("request_id", "unknown")
            print(f"[agent] Cancel: {request_id}")
            await self.interrupt(msg)
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
        elif msg_type == "answer":
            print(f"[agent] {msg.get('answers')}")

        elif msg_type == "kernel_message":
            nested_msg = msg.get("message", {})
            nested_type = nested_msg.get("type")
            print(f"[agent] Kernel message: {nested_type}")

            if nested_type == "cell_result":
                cell_id = nested_msg.get("cell_id")
                success = not nested_msg.get("has_exception", False)
                cell_name = nested_msg.get("display_name")

                logs = nested_msg.get("logs", None)
                if logs is not None and len(logs) > 4096:
                    logs = logs[-4096:]

                if cell_id not in self.pending_cells:
                    print(f"[agent] Ignoring cell_result for {cell_id} ")
                    return

                if cell_id is not None:
                    self.pending_cells.discard(str(cell_id))

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
                    exception = nested_msg.get("exception", "Unknown error")
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
                assert self.claude is not None

                prompt_content = f"Cell execution update:\n\n {result_message}"
                if not success:
                    exception = result_content.get("exception")
                    if exception:
                        prompt_content = f"{prompt_content}\n\nException: {exception}"
                if logs:
                    prompt_content = f"{prompt_content}\n\nLogs:\n{logs}"

                await self.query(
                    AgentQuery(
                        type="agent_query",
                        query=prompt_content,
                        behavior=self.behavior,
                        template_version_id=None,
                    )
                )

            elif nested_type == "set_widget_value":
                data = nested_msg.get("data", {})
                for key, value in data.items():
                    if key in self.pending_widgets:
                        self.pending_widgets[key] = value
                        print(f"        Set widget {key}")

                if all(v is not None for v in self.pending_widgets.values()):
                    await self.set_agent_status("thinking")

                    widget_info = ", ".join(
                        f"{k}={v}" for k, v in self.pending_widgets.items()
                    )
                    content = f"User provided input via widget(s): {widget_info}"

                    await self._insert_history(
                        payload={"content": content, "hidden": True}
                    )
                    await self.query(
                        AgentQuery(
                            type="agent_query",
                            query=content,
                            behavior=self.behavior,
                            template_version_id=None,
                        )
                    )
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
                await harness.disconnect()
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
