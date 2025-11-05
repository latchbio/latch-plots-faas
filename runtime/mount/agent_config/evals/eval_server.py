import argparse
import asyncio
import json
import os
import socket
import sys
import textwrap
import time
from pathlib import Path

os.environ.setdefault("auto_reload", "false")
os.environ.setdefault("logging_mode", "console")
os.environ.setdefault("domain", "latch.bio")
os.environ.setdefault("DD_VERSION", "eval")
os.environ.setdefault("DD_SERVICE", "latch-plots-eval")
os.environ.setdefault("DD_ENV", "eval")
os.environ.setdefault("DD_AGENT_HOST", "localhost")
os.environ.setdefault("DD_TRACE_ENABLED", "false")
os.environ.setdefault("DD_PROFILING_ENABLED", "false")
os.environ.setdefault("DD_RUNTIME_METRICS_ENABLED", "false")
os.environ.setdefault("OTEL_SDK_DISABLED", "true")

import websockets
from eval_types import TestCase, TestResult
from run_local_eval import get_eval_config
from binary_grader import GRADER_REGISTRY

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent))

from mount.socketio import SocketIo
from mount.utils import gql_query


class EvalServer:
    def __init__(self, test_case: TestCase, latch_dir: Path):
        self.test_case = test_case
        self.latch_dir = latch_dir
        self.conversation_history = []
        self.last_message_time = None
        self.agent_sent_result = False
        self.has_questions = False
        self.cell_status = {}
        self.agent_proc = None
        self.agent_sock = None
        self.agent_conn = None
        self.websocket = None
        self.notebook_context = None
        self.cleanup_complete = False
        self.session_id = None
        self.eval_complete = False

    async def start_agent(self):
        print(f"[eval] Starting agent for test: {self.test_case.id}")

        self.session_id = None

        sock_a, sock_agent = socket.socketpair(family=socket.AF_UNIX)
        sock_a.setblocking(False)
        sock_agent_fd = sock_agent.detach()

        self.agent_sock = sock_a
        self.agent_conn = await SocketIo.from_socket(sock_a)

        agent_path = Path(__file__).parent.parent.parent / "agent.py"

        self.agent_proc = await asyncio.create_subprocess_exec(
            sys.executable,
            "-u",
            str(agent_path),
            str(sock_agent_fd),
            pass_fds=[sock_agent_fd],
            stdin=asyncio.subprocess.DEVNULL,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            env={
                **os.environ,
                "LATCH_SANDBOX_ROOT": str(self.latch_dir),
                "PYTHONUNBUFFERED": "1",
                "AGENT_DEBUG": "1",
            },
            preexec_fn=lambda: os.nice(5),
            limit=1024 * 1024,
        )

        async def stream_output(stream, prefix=""):
            while True:
                try:
                    line = await stream.readline()
                    if not line:
                        break
                    decoded = line.decode().rstrip()
                    if len(decoded) > 1000:
                        decoded = decoded[:1000] + "... [TRUNCATED]"
                    print(f"[agent] {prefix}{decoded}", flush=True)
                except (ValueError, asyncio.LimitOverrunError) as e:
                    if "limit" in str(e).lower():
                        chunk = await stream.read(8192)
                        if not chunk:
                            break
                        print(f"[agent] {prefix}[Large output truncated: {len(chunk)} bytes]", flush=True)
                    else:
                        raise
                except Exception as e:
                    print(f"[agent] {prefix}[Error reading output: {e}]", flush=True)
                    break

        asyncio.create_task(stream_output(self.agent_proc.stdout, ""))
        asyncio.create_task(stream_output(self.agent_proc.stderr, "[stderr] "))

        msg = await self.agent_conn.recv()
        if msg.get("type") == "ready":
            print("[eval] Agent subprocess started and ready")

    async def _send_action_request(self, action: str, params: dict = {}) -> dict | None:
        if not self.websocket:
            print(f"[eval] Error: websocket is None for action {action}")
            return None

        tx_id = f"eval_{action}_{time.time()}"
        request = {
            "type": "agent_action",
            "action": action,
            "params": params,
            "tx_id": tx_id
        }

        try:
            await self.websocket.send(json.dumps(request))
        except Exception as e:
            print(f"[eval] Error sending {action} request: {e}")
            return None

        try:
            response_msg = await asyncio.wait_for(self.websocket.recv(), timeout=5.0)
            response = json.loads(response_msg)

            if response.get("tx_id") == tx_id and response.get("type") == "agent_action_response":
                if response.get("status") == "success":
                    return response
                else:
                    print(f"[eval] Failed {action}: status={response.get('status')}, error={response.get('error')}")
                    return None
            else:
                print(f"[eval] Unexpected response type: {response.get('type')}")
                return None
        except asyncio.TimeoutError:
            print(f"[eval] Timeout waiting for {action}")
            return None
        except Exception as e:
            print(f"[eval] Error getting {action}: {e}")
            return None

    async def get_notebook_context_from_console(self):
        print("[eval] Requesting notebook context from console...")

        context_result = await self._send_action_request("get_context")

        result = {}

        if isinstance(context_result, dict) and context_result is not None:
            result["context"] = context_result.get("context", {})
        else:
            print(f"[eval] Failed to get notebook context: {context_result}")
            result["context"] = {}

        print(f"[eval] Got notebook context: {len(result.get('context', {}).get('cells', []))} cells")

        return result

    async def stop_agent(self):
        if self.agent_proc:
            print("[eval] Stopping agent")
            try:
                self.agent_proc.terminate()
                await asyncio.wait_for(self.agent_proc.wait(), timeout=2)
            except TimeoutError:
                self.agent_proc.kill()
                await self.agent_proc.wait()

        if self.agent_sock:
            try:
                self.agent_sock.close()
            except Exception:
                pass

        self.agent_proc = None
        self.agent_sock = None
        self.agent_conn = None

    async def handle_agent_connection(self, websocket):
        print("[eval] Agent websocket connected")
        self.websocket = websocket

        try:
            # Wait for console to send init with its session_id, then use that
            print("[eval] Waiting for console init to get session_id...")
            init_msg = await websocket.recv()
            console_init = json.loads(init_msg)
            if console_init.get("type") == "init":
                self.session_id = int(console_init.get("session_id"))
    
                # Clear history for this session
                auth_token_sdk, _, _ = get_eval_config()

                try:
                    await gql_query(
                        auth=auth_token_sdk,
                        query="""
                            mutation ClearAgentHistory($sessionId: BigInt!) {
                                clearAgentHistory(input: {argSessionId: $sessionId}) {
                                    clientMutationId
                                }
                            }
                        """,
                        variables={"sessionId": str(self.session_id)},
                    )
                except Exception as e:
                    print(f"[eval] Warning: Failed to clear agent history: {e}")
            await self.agent_conn.send({"type": "init", "session_id": self.session_id})

            msg = await self.agent_conn.recv()
            if msg.get("type") == "agent_status" and msg.get("status") == "ready":
                data_context = ""
                if self.test_case.data_node:
                    data_nodes = self.test_case.data_node if isinstance(self.test_case.data_node, list) else [self.test_case.data_node]
                    contextual_data = []
                    for node in data_nodes:
                        contextual_data.append({
                            "type": "File",
                            "path": node,
                            "id": node.replace("latch:///", "").replace(".csv", "").replace(".h5ad", ""),
                        })
                    data_context = f"\n\nHere is the context of the selected nodes the user would like to use: <ContextualNodeData>{json.dumps(contextual_data)}</ContextualNodeData>"

                initial_query = textwrap.dedent(f"""
                    First, delete all cells in the notebook to start fresh. Then: {self.test_case.task}

                    IMPORTANT: When you have completed this task, include the text "[EVAL_COMPLETE]" at the end of your summary in submit_response.
                    {data_context}
                """).strip()

                await self.agent_conn.send({
                    "type": "agent_query",
                    "query": initial_query,
                    "request_id": f"eval-init-{self.session_id}"
                })

            async def forward_agent_to_console():
                try:
                    while True:
                        msg = await self.agent_conn.recv()
                        msg_type = msg.get("type", "unknown")
                        print(f"[eval] agent→console: {msg_type}")
                        await self.handle_agent_message(msg)
                        await websocket.send(json.dumps(msg))
                except Exception as e:
                    print(f"[eval] Error forwarding agent→console: {e}")
                    import traceback
                    traceback.print_exc()

            async def forward_console_to_agent():
                try:
                    async for message in websocket:
                        msg = json.loads(message)
                        msg_type = msg.get("type", "unknown")
                        print(f"[eval] console→agent: {msg_type}")
                        await self.agent_conn.send(msg)
                except Exception as e:
                    print(f"[eval] Error forwarding console→agent: {e}")
                    import traceback
                    traceback.print_exc()

            forward_task = asyncio.create_task(forward_agent_to_console())
            receive_task = asyncio.create_task(forward_console_to_agent())

            while not self.is_done():
                if forward_task.done() or receive_task.done():
                    print("[eval] One of the forwarding tasks completed unexpectedly")
                    if forward_task.done():
                        try:
                            forward_task.result()
                        except Exception as e:
                            print(f"[eval] Forward task error: {e}")
                    if receive_task.done():
                        try:
                            receive_task.result()
                        except Exception as e:
                            print(f"[eval] Receive task error: {e}")
                    break
                await asyncio.sleep(1)

            print("[eval] Test complete, stopping console→agent forwarding...")
            receive_task.cancel()
            try:
                await asyncio.wait_for(receive_task, timeout=0.1)
            except (asyncio.CancelledError, asyncio.TimeoutError):
                pass

            print("[eval] Getting notebook context...")
            self.notebook_context = await self.get_notebook_context_from_console()

            if not self.notebook_context:
                print("[eval] Warning: Failed to get notebook context, websocket may be closed")

            print("[eval] Stopping agent→console forwarding...")
            forward_task.cancel()
            try:
                await asyncio.wait_for(forward_task, timeout=0.1)
            except (asyncio.CancelledError, asyncio.TimeoutError):
                pass

            await self.stop_agent()

            print("[eval] Cleanup complete")
            self.cleanup_complete = True
        except Exception as e:
            print(f"[eval] Error in handle_agent_connection: {e}")
            import traceback
            traceback.print_exc()
            await self.stop_agent()
            self.cleanup_complete = True
            raise

    async def handle_agent_message(self, msg: dict):
        self.conversation_history.append(msg)
        self.last_message_time = time.time()

        msg_type = msg.get("type")

        if msg_type == "agent_history_updated":
            await self.check_for_completion()

        elif msg_type in ("start_cell", "cell_result"):
            cell_id = msg.get("cell_id")
            if cell_id:
                if msg_type == "start_cell":
                    self.cell_status[cell_id] = "running"
                elif msg_type == "cell_result":
                    has_exception = msg.get("has_exception", False)
                    self.cell_status[cell_id] = "error" if has_exception else "ran"

    async def check_for_completion(self):
        auth_token_sdk, _, _ = get_eval_config()

        try:
            resp = await gql_query(
                auth=auth_token_sdk,
                query="""
                    query AgentHistory($sessionId: BigInt!) {
                        agentHistories(condition: {sessionId: $sessionId, removed: false}, orderBy: ID_DESC, first: 20) {
                            nodes { id payload }
                        }
                    }
                """,
                variables={"sessionId": str(self.session_id)},
            )
            nodes = resp.get("data", {}).get("agentHistories", {}).get("nodes", [])

            for node in nodes:
                payload = node.get("payload", {})
                if payload.get("type") == "anthropic_message" and payload.get("role") == "assistant":
                    content = payload.get("content", [])
                    for block in content:
                        if isinstance(block, dict) and block.get("type") == "tool_use" and block.get("name") == "submit_response":
                            tool_input = block.get("input", {})
                            summary = tool_input.get("summary", "")
                            if "[EVAL_COMPLETE]" in summary:
                                self.eval_complete = True
                                return
        except Exception as e:
            print(f"[eval] Error checking for completion: {e}")

    async def fetch_full_conversation_history(self) -> list[dict]:
        auth_token_sdk, _, _ = get_eval_config()

        try:
            resp = await gql_query(
                auth=auth_token_sdk,
                query="""
                    query AgentHistory($sessionId: BigInt!) {
                        agentHistories(condition: {sessionId: $sessionId, removed: false}, orderBy: ID_ASC) {
                            nodes { id payload }
                        }
                    }
                """,
                variables={"sessionId": str(self.session_id)},
            )
            nodes = resp.get("data", {}).get("agentHistories", {}).get("nodes", [])
            return [node.get("payload", {}) for node in nodes]
        except Exception as e:
            print(f"[eval] Error fetching conversation history: {e}")
            return []

    def is_done(self) -> bool:
        if not self.eval_complete:
            return False

        if self.has_questions:
            return False

        running_cells = [cid for cid, status in self.cell_status.items() if status == "running"]
        if running_cells:
            print(f"[eval] Not done: {len(running_cells)} cells still running: {running_cells}")
            return False

        print(f"[eval] Completion detected: [EVAL_COMPLETE] marker found, no questions, no running cells")
        return True


async def run_eval(test_case: TestCase, port: int, latch_dir: Path) -> TestResult:
    print(f"\n{'=' * 70}")
    print(f"Running eval: {test_case.id}")
    print(f"Listening on port: {port}")
    print("=" * 70)

    start_time = time.time()

    server = EvalServer(test_case, latch_dir)

    async def connection_handler(websocket, path):
        if path == "/agent":
            print("[eval] Console connected")
            await server.start_agent()
            await server.handle_agent_connection(websocket)
        else:
            print(f"[eval] Unknown path: {path}")
            await websocket.close()

    async with websockets.serve(
        connection_handler,
        "localhost",
        port,
        max_size=10 * 1024 * 1024
    ):
        print(f"[eval] WebSocket server listening on ws://localhost:{port}/agent")
        print("[eval] Waiting for console to connect to agent endpoint...")
        print("[eval] Note: Console should connect to remote kernel separately")

        while not server.cleanup_complete:
            await asyncio.sleep(1)

    duration_ms = (time.time() - start_time) * 1000

    if not server.notebook_context:
        error_msg = "Failed to retrieve notebook context from console. "
        if not server.websocket:
            error_msg += "Console never connected to the eval server. "
        error_msg += "Make sure the console is connected to ws://localhost:{port}/agent before starting the eval."
        raise RuntimeError(error_msg)

    notebook_state = server.notebook_context

    print("[eval] Fetching full conversation history from database...")
    conversation_history = await server.fetch_full_conversation_history()

    test_result = TestResult(
        test_id=test_case.id,
        conversation_history=conversation_history,
        notebook_state=notebook_state,
        duration_ms=duration_ms,
    )

    if test_case.grader:
        print("[eval] Running binary grader...")
        grader_type = test_case.grader.get("type")
        grader_config = test_case.grader.get("config", {})

        if grader_type in GRADER_REGISTRY:
            grader_cls = GRADER_REGISTRY[grader_type]
            grader = grader_cls()
            grader_result = grader.evaluate(test_result, grader_config)

            test_result.grader_result = {
                "passed": grader_result.passed,
                "metrics": grader_result.metrics,
                "reasoning": grader_result.reasoning,
                "agent_answer": grader_result.agent_answer
            }

            print(f"[eval] Grader result: {'PASS' if grader_result.passed else 'FAIL'}")
            print(f"[eval] Grader reasoning:\n{grader_result.reasoning}")
        else:
            print(f"[eval] Warning: Unknown grader type '{grader_type}'")

    print(f"\n[eval] Eval completed in {duration_ms / 1000:.2f}s")
    print(f"[eval] Total conversation turns: {len(conversation_history)}")
    cells = notebook_state.get("context", {}).get("cells", [])
    print(f"[eval] Final notebook state: {len(cells)} cells")

    return test_result


async def main():
    parser = argparse.ArgumentParser(description="Run agent eval server")
    parser.add_argument("--eval", required=True, help="Eval file to run")
    parser.add_argument("--port", type=int, default=8765, help="WebSocket server port")
    args = parser.parse_args()

    eval_dir = Path(__file__).parent
    eval_file = eval_dir / args.eval

    with open(eval_file) as f:
        test_data = json.load(f)

    test_case = TestCase(**test_data)

    import shutil
    latch_dir = Path.home() / ".latch" / "eval" / test_case.id

    if latch_dir.exists():
        print(f"[eval] Removing existing sandbox at {latch_dir}")
        shutil.rmtree(latch_dir)

    latch_dir.mkdir(parents=True, exist_ok=True)
    print(f"[eval] Created fresh sandbox at {latch_dir}")

    user_token = Path.home() / ".latch" / "token"
    if user_token.exists():
        (latch_dir / "token").write_text(user_token.read_text())
    else:
        (latch_dir / "token").write_text("local-dev-token")

    (latch_dir / "id").write_text("99999")
    (latch_dir / "session-id").write_text(f"eval-{test_case.id}")
    (latch_dir / "nucleus-url").write_text("https://nucleus.latch.bio")

    result = await run_eval(test_case, args.port, latch_dir)

    output_file = eval_dir / f"result_{test_case.id}.json"
    with open(output_file, "w") as f:
        json.dump(result.model_dump(), f, indent=2)
    print(f"\n[eval] Result saved to {output_file}")


if __name__ == "__main__":
    asyncio.run(main())
