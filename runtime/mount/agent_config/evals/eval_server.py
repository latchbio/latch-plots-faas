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
os.environ.setdefault("domain", "localhost")
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
from judge import LLMJudge
from run_local_eval import get_eval_config

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent))

from mount.socketio import SocketIo


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
        self.kernel_proc = None
        self.kernel_sock = None
        self.kernel_conn = None
        self.websocket = None
        self.kernel_websocket = None
        self.notebook_context = None
        self.cleanup_complete = False
        self.session_id = None
        self.eval_complete = False

    async def start_kernel(self):
        print("[eval] Starting kernel subprocess")

        sock_k, sock_kernel = socket.socketpair(family=socket.AF_UNIX)
        sock_k.setblocking(False)
        sock_kernel_fd = sock_kernel.detach()

        self.kernel_sock = sock_k
        self.kernel_conn = await SocketIo.from_socket(sock_k)

        kernel_path = Path(__file__).parent.parent.parent / "kernel.py"

        self.kernel_proc = await asyncio.create_subprocess_exec(
            sys.executable,
            str(kernel_path),
            str(sock_kernel_fd),
            pass_fds=[sock_kernel_fd],
            stdin=asyncio.subprocess.DEVNULL,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            env={
                **os.environ,
                "LATCH_SANDBOX_ROOT": str(self.latch_dir),
                "PYTHONUNBUFFERED": "1",
            },
            preexec_fn=lambda: os.nice(1),
        )

        async def stream_output(stream, prefix=""):
            while True:
                line = await stream.readline()
                if not line:
                    break
                print(f"[kernel] {prefix}{line.decode().rstrip()}", flush=True)

        asyncio.create_task(stream_output(self.kernel_proc.stdout, ""))
        asyncio.create_task(stream_output(self.kernel_proc.stderr, "[stderr] "))

        msg = await self.kernel_conn.recv()
        if msg.get("type") == "ready":
            print("[eval] Kernel ready, sending init...")
            await self.kernel_conn.send({
                "type": "init",
                "widget_states": {},
                "cell_output_selections": {},
                "plot_data_selections": {},
                "viewer_cell_data": {},
                "plot_configs": {},
                "session_snapshot_mode": False,
            })
            print("[eval] Kernel initialized")

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
        )

        async def stream_output(stream, prefix=""):
            while True:
                line = await stream.readline()
                if not line:
                    break
                print(f"[agent] {prefix}{line.decode().rstrip()}", flush=True)

        asyncio.create_task(stream_output(self.agent_proc.stdout, ""))
        asyncio.create_task(stream_output(self.agent_proc.stderr, "[stderr] "))

        msg = await self.agent_conn.recv()
        if msg.get("type") == "ready":
            print("[eval] Agent subprocess started and ready")

    async def get_notebook_context_from_console(self):
        print("[eval] Requesting notebook context from console...")

        if not self.websocket:
            print("[eval] Error: websocket is None")
            return None

        tx_id = "eval_get_context"
        context_request = {
            "type": "agent_action",
            "action": "get_context",
            "params": {},
            "tx_id": tx_id
        }

        try:
            await self.websocket.send(json.dumps(context_request))
            print("[eval] Sent context request, waiting for response...")
        except Exception as e:
            print(f"[eval] Error sending context request: {e}")
            return None

        try:
            response_msg = await asyncio.wait_for(self.websocket.recv(), timeout=5.0)
            response = json.loads(response_msg)

            if response.get("tx_id") == tx_id and response.get("type") == "agent_action_response":
                if response.get("status") == "success":
                    context = response.get("context", {})
                    return context
                else:
                    print(f"[eval] Failed to get context: status={response.get('status')}, error={response.get('error')}")
                    return None
            else:
                print(f"[eval] Unexpected response type: {response.get('type')}, expected 'agent_action_response'")
                return None
        except asyncio.TimeoutError:
            print("[eval] Timeout waiting for notebook context")
            return None
        except Exception as e:
            print(f"[eval] Error getting notebook context: {e}")
            import traceback
            traceback.print_exc()
            return None

    async def stop_kernel(self):
        if self.kernel_proc:
            print("[eval] Stopping kernel")
            try:
                self.kernel_proc.terminate()
                await asyncio.wait_for(self.kernel_proc.wait(), timeout=2)
            except TimeoutError:
                self.kernel_proc.kill()
                await self.kernel_proc.wait()

        if self.kernel_sock:
            try:
                self.kernel_sock.close()
            except Exception:
                pass

        self.kernel_proc = None
        self.kernel_sock = None
        self.kernel_conn = None

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

    async def handle_kernel_connection(self, websocket):
        print("[eval] Kernel websocket connected")
        self.kernel_websocket = websocket

        try:
            auth_msg = await websocket.recv()
            auth_data = json.loads(auth_msg)
            print(f"[eval] Received auth message: notebook_id={auth_data.get('notebook_id')}, token={auth_data.get('token')[:20]}...")

            print("[eval] Sending ready message to console...")
            await websocket.send(json.dumps({
                "type": "ready",
                "connection_idx": 0,
                "cell_status": {},
                "cell_sequencers": {},
                "cell_outputs": {}
            }))

            async def forward_kernel_to_console():
                try:
                    while not self.cleanup_complete:
                        msg = await self.kernel_conn.recv()
                        await websocket.send(json.dumps(msg))
                except Exception as e:
                    if not self.cleanup_complete:
                        print(f"[eval] Error forwarding kernel→console: {e}")
                        import traceback
                        traceback.print_exc()

            async def forward_console_to_kernel():
                try:
                    async for message in websocket:
                        if self.cleanup_complete:
                            break
                        msg = json.loads(message)
                        await self.kernel_conn.send(msg)
                except Exception as e:
                    if not self.cleanup_complete:
                        print(f"[eval] Error forwarding console→kernel: {e}")
                        import traceback
                        traceback.print_exc()

            forward_task = asyncio.create_task(forward_kernel_to_console())
            receive_task = asyncio.create_task(forward_console_to_kernel())

            while not self.cleanup_complete:
                await asyncio.sleep(0.1)

            forward_task.cancel()
            receive_task.cancel()
            try:
                await asyncio.gather(forward_task, receive_task, return_exceptions=True)
            except:
                pass

            try:
                await websocket.close()
            except:
                pass
        except Exception as e:
            print(f"[eval] Error in handle_kernel_connection: {e}")
            import traceback
            traceback.print_exc()

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
                from run_local_eval import get_eval_config
                from mount.utils import gql_query
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
                    contextual_data = [{
                        "type": "File",
                        "path": self.test_case.data_node,
                        "id": self.test_case.data_node.replace("latch:///", "").replace(".csv", ""),
                    }]
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
            await self.stop_kernel()

            print("[eval] Cleanup complete")
            self.cleanup_complete = True
        except Exception as e:
            print(f"[eval] Error in handle_agent_connection: {e}")
            import traceback
            traceback.print_exc()
            await self.stop_agent()
            await self.stop_kernel()
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
        from run_local_eval import get_eval_config
        from mount.utils import gql_query
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

    kernel_connected = asyncio.Event()
    agent_connected = asyncio.Event()

    async def connection_handler(websocket, path):
        if path == "/run":
            print("[eval] Kernel connection received, starting kernel subprocess...")
            await server.start_kernel()
            kernel_connected.set()
            await server.handle_kernel_connection(websocket)
        elif path == "/agent":
            print("[eval] Agent connection received, waiting for kernel...")
            await kernel_connected.wait()
            print("[eval] Starting agent subprocess...")
            await server.start_agent()
            agent_connected.set()
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
        print(f"[eval] WebSocket server listening on ws://localhost:{port}")
        print(f"[eval]   Kernel endpoint: ws://localhost:{port}/run")
        print(f"[eval]   Agent endpoint: ws://localhost:{port}/agent")
        print("[eval] Waiting for console to connect...")

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

    test_result = TestResult(
        test_id=test_case.id,
        conversation_history=server.conversation_history,
        notebook_state=notebook_state,
        duration_ms=duration_ms,
    )

    print(f"\n[eval] Eval completed in {duration_ms / 1000:.2f}s")
    print(f"[eval] Total conversation turns: {len(server.conversation_history)}")
    cells = notebook_state.get("cells", [])
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

    auth_token_sdk, nucleus_url, pod_id = get_eval_config()
    judge = LLMJudge(
        api_key="dummy",
        base_url=f"{nucleus_url}/infer/plots-agent/anthropic",
        headers={"Authorization": auth_token_sdk, "Pod-Id": str(pod_id)}
    )

    print("\n[eval] Judging result...")
    eval_result = await judge.evaluate(test_case, result)
    result.eval_result = eval_result

    print(f"\n{'=' * 70}")
    print(f"Score: {eval_result.score:.2f}")
    print(f"Passed: {'✓ PASS' if eval_result.passed else '✗ FAIL'}")
    print("\nReasoning:")
    print(eval_result.reasoning)
    print(f"\nSuccesses ({len(eval_result.successes)}):")
    for success in eval_result.successes:
        print(f"  ✓ {success}")
    print(f"\nFailures ({len(eval_result.failures)}):")
    for failure in eval_result.failures:
        print(f"  ✗ {failure}")
    print("=" * 70)

    output_file = eval_dir / f"result_{test_case.id}.json"
    with open(output_file, "w") as f:
        json.dump(result.model_dump(), f, indent=2)
    print(f"\n[eval] Result saved to {output_file}")


if __name__ == "__main__":
    asyncio.run(main())
