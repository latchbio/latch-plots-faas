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
    HookContext,
    McpSdkServerConfig,
    PostToolUseHookInput,
    StreamEvent,
    SyncHookJSONOutput,
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


Behavior = Literal["plan", "build"]

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


ToolResultContent = str | list[TextToolResultBlock]


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


class AgentQuery(TypedDict):
    type: Literal["agent_query"]
    request_id: str
    query: str  # todo(rteqs): rename to prompt
    behavior: Behavior
    template_version_id: str
    contextual_node_data: NotRequired[dict[str, Any]]  # todo(rteqs): type properly
    selected_widgets: NotRequired[list[dict[str, Any]]]  # todo(rteqs): type properly


@dataclass
class AgentHarness:
    conn: SocketIoThread
    claude: ClaudeSDKClient | None = None
    system_prompy: str | None = None
    pending_operations: dict[str, asyncio.Future] = field(default_factory=dict)
    executing_cells: set[str] = field(default_factory=set)
    operation_counter: int = 0
    current_request_id: str | None = None

    conversation_task: asyncio.Task | None = None
    agent_session_id: int | None = None  # db_session_id
    claude_session_id: str | None = None
    latest_notebook_context: dict = field(default_factory=dict)
    current_status: str | None = None
    expected_widgets: dict[str, Any | None] = field(default_factory=dict)
    behavior: Behavior = "plan"
    latest_notebook_state: str | None = None
    current_plan: dict | None = None
    in_memory_history: list[dict] = field(default_factory=list)
    mcp_server: McpSdkServerConfig = field(default_factory=lambda: agent_tools_mcp)
    mcp_allowed_tools: list[str] = field(
        default_factory=lambda: list(MCP_ALLOWED_TOOL_NAMES)
    )
    current_query_task: asyncio.Task | None = None

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

    # todo(rteqs): for old sessions that need to be migrated. we need to inject out session history from vacuole into the .jsonl history file

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

    async def _load_claude_session_id(self, session_id: int) -> str | None:
        # todo(rteqs): still cringe
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

    # todo(rteqs): make reliable / fault tolerant
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
                # todo(rteqs): claude code + our system prompt
                system_prompt=self.system_prompt,
                include_partial_messages=False,
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
                can_use_tool=can_use_tool,
            )
        )

        tries = 3
        for _ in range(tries):
            try:
                await self.claude.connect()
                break

                # server_info = await self.claude.get_server_info()
                # print(server_info)
                # if server_info is None:
                #     continue

                # session_id = server_info.get("session_id")
                # if session_id is not None:
                #     self.claude_session_id = session_id
                # break

            except Exception:
                traceback.print_exc()
        else:
            print("[agent] Failed to connect ")
            raise RuntimeError("failed to connect")

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
            await self.disconnect()

        resume_session_id = msg.get("claude_session_id")
        if resume_session_id is None or resume_session_id == "":
            resume_session_id = await self._load_claude_session_id(
                self.agent_session_id
            )
        self.claude_session_id = resume_session_id

        if self.claude is None:
            await self.connect(resume_session_id=resume_session_id)
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

    async def _reset_for_new_session(self) -> None:
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

        self.executing_cells.clear()
        self.expected_widgets.clear()
        self.current_request_id = None
        self.current_status = None
        self.current_plan = None

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
            # todo(rteqs): send to browser to start stream?
            return

        if event.type == "message_delta":
            # todo(rteqs): handle stop_reason
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

    async def create_prompt(self, query: AgentQuery) -> str:
        # todo(rteqs): implement
        # logic from _wait_for_message, _run_query_with_trun_prompt, _run_query here
        return query["query"]

    async def query(self, msg: AgentQuery) -> None:
        # todo(rteqs): implement logic for wait for connnetion
        assert self.claude is not None
        # assert self.claude_session_id is not None

        prompt = await self.create_prompt(msg)

        await self.send({
            "type": "agent_stream_start",
            "timestamp": int(time.time() * 1000),
        })

        print(
            "[agent] starting SDK query "
            f"(request_id={msg['request_id']}, db_session_id={self.agent_session_id}, "
            f"claude_session_id={self.claude_session_id})"
        )

        await self._insert_history(
            role="user",
            request_id=msg["request_id"],
            payload={"content": msg["query"], "display_query": msg["query"]},
        )
        await self.claude.query(prompt=prompt, session_id=self.claude_session_id)

        async for res in self.claude.receive_response():
            # todo(rteqs): pretend we can't be interrupted for now
            if isinstance(res, StreamEvent):
                await self._handle_stream_event(res)

            elif isinstance(res, AssistantMessage):
                if res.error is not None:
                    # todo(rteqs): do we care here? doesn't ResultMessage give us the error also?
                    print(f"[agent] assistant message error={res.error}")

                for c in res.content:
                    if isinstance(c, TextBlock):
                        await self._insert_history(
                            role="assistant",
                            request_id=msg["request_id"],
                            payload={"content": [{"type": "text", "text": c.text}]},
                        )

                    elif isinstance(c, ThinkingBlock):
                        await self._insert_history(
                            role="assistant",
                            request_id=msg["request_id"],
                            payload={
                                "content": [
                                    {"type": "thinking", "thinking": c.thinking}
                                ]
                            },
                        )

                    elif isinstance(c, ToolUseBlock):
                        await self._insert_history(
                            role="assistant",
                            request_id=msg["request_id"],
                            payload={
                                "content": [
                                    {
                                        "type": "tool_use",
                                        "id": c.id,
                                        "name": c.name,
                                        "input": c.input,
                                    }
                                ]
                            },
                        )

            elif isinstance(res, UserMessage):
                if isinstance(res.content, str):
                    await self._insert_history(
                        role="user",
                        request_id=msg["request_id"],
                        payload={"content": res.content},
                    )
                    continue

                for c in res.content:
                    if isinstance(c, TextBlock):
                        await self._insert_history(
                            role="user",
                            request_id=msg["request_id"],
                            payload={"content": [{"type": "text", "text": c.text}]},
                        )

                    if isinstance(c, ToolResultBlock):
                        # todo(rteqs): update console AgentToolResult to have is_error
                        # todo(rteqs): pass c.content to get GenericToolResult
                        await self._insert_history(
                            role="user",
                            request_id=msg["request_id"],
                            payload={
                                "content": [
                                    {
                                        "type": "tool_result",
                                        "tool_use_id": c.tool_use_id,
                                        "content": c.content,
                                        "is_error": c.is_error,
                                    }
                                ]
                            },
                        )

            elif isinstance(res, ResultMessage):
                # if res.session_id != self.claude_session_id:
                #     continue

                usage = validate(res.usage, Usage)
                await self._send_usage_update(usage)

                # todo(rteqs): do we need to send a message frontend or history to say we're done?
                # if res.
                if res.subtype == "success":
                    break

                assert res.result is None
                assert res.is_error is True

                # todo(rteqs): self.handle_error_result

                break

            elif isinstance(res, SystemMessage):
                ...
                # todo(rteqs): handle this for processing subagent messages
                # note: we should just ignore init messages. nothing important and we shouldn't be changing our session_id mid query

    async def interrupt(self, msg: dict[str, Any]) -> None:
        # todo(rteqs): clean up
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

        messages = None
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
        self.system_prompt = new_content
        self.claude.options.system_prompt = new_content

        full_prompt = await self.get_full_prompt()
        return {"status": "success", **full_prompt}

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
            query_preview = query[:60] + "…" if len(query) > 60 else query
            request_id = msg.get("request_id", "unknown")
            print(
                f"[agent] accept: dispatching to handle_query (query={query_preview}, request_id={request_id})"
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
        elif msg_type == "kernel_message":
            nested_msg = msg.get("message", {})
            nested_type = nested_msg.get("type")
            print(f"[agent] Kernel message: {nested_type}")

            # todo(rteqs): fill in logic from wait_for_message for cell_result and set_widget_value
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

                # if self.current_request_id is not None:
                #     await self.pending_messages.put({
                #         "type": "cell_result",
                #         "cell_id": cell_id,
                #         "success": not has_exception,
                #         "exception": exception,
                #         "logs": logs,
                #         "display_name": display_name,
                #     })
                # else:
                #     print(
                #         f"        Cell {cell_id} completed but no active request - updating executing_cells only"
                #     )

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
                        # await self.pending_messages.put({
                        #     "type": "set_widget_value",
                        #     "data": self.expected_widgets,
                        # })
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
