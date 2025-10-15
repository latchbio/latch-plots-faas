import argparse
import asyncio
import json
import os
import socket
import sys
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
        self.websocket = None

    async def start_agent(self):
        print(f"[eval] Starting agent for test: {self.test_case.id}")

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
            "--initial-query",
            self.test_case.task,
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

        await self.agent_conn.send({"type": "init"})

        msg = await self.agent_conn.recv()
        if msg.get("type") == "ready":
            print("[eval] Agent ready")

    async def stop_agent(self):
        if self.agent_proc:
            print("[eval] Stopping agent")
            try:
                self.agent_proc.terminate()
                await asyncio.wait_for(self.agent_proc.wait(), timeout=2)
            except asyncio.TimeoutError:
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

    async def handle_console_connection(self, websocket):
        print("[eval] Console connected")
        self.websocket = websocket

        try:
            await self.start_agent()

            async def forward_agent_to_console():
                try:
                    while True:
                        msg = await self.agent_conn.recv()
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

            print("[eval] Test complete, stopping agent")
            forward_task.cancel()
            receive_task.cancel()

            await self.stop_agent()
        except Exception as e:
            print(f"[eval] Error in handle_console_connection: {e}")
            import traceback
            traceback.print_exc()
            await self.stop_agent()
            raise

    async def handle_agent_message(self, msg: dict):
        self.conversation_history.append(msg)
        self.last_message_time = time.time()

        msg_type = msg.get("type")

        if msg_type == "agent_result":
            self.agent_sent_result = True
            structured = msg.get("structured_output", {})
            questions = structured.get("questions")
            self.has_questions = questions is not None and len(questions) > 0

            if self.has_questions:
                print(f"[eval] Agent asked {len(questions)} question(s), auto-responding...")
                await asyncio.sleep(1)
                auto_response = {
                    "type": "agent_query",
                    "query": "Please continue and make your best judgment.",
                    "request_id": f"auto-response-{time.time()}"
                }
                await self.agent_conn.send(auto_response)
                self.has_questions = False
                self.agent_sent_result = False
                print(f"[eval] Waiting for agent to respond to auto-answer...")
            else:
                print(f"[eval] Agent sent result (no questions)")

        elif msg_type in ("start_cell", "cell_result"):
            cell_id = msg.get("cell_id")
            if cell_id:
                if msg_type == "start_cell":
                    self.cell_status[cell_id] = "running"
                elif msg_type == "cell_result":
                    has_exception = msg.get("has_exception", False)
                    self.cell_status[cell_id] = "error" if has_exception else "ran"

    def is_done(self) -> bool:
        if not self.agent_sent_result:
            return False

        if self.has_questions:
            return False

        running_cells = [cid for cid, status in self.cell_status.items() if status == "running"]
        if running_cells:
            return False

        if self.last_message_time is None:
            return False

        idle_time = time.time() - self.last_message_time
        if idle_time < 10:
            return False

        print(f"[eval] Completion detected: result sent, no questions, no running cells, idle for {idle_time:.1f}s")
        return True

async def run_eval(test_case: TestCase, port: int, latch_dir: Path) -> TestResult:
    print(f"\n{'=' * 70}")
    print(f"Running eval: {test_case.id}")
    print(f"Listening on port: {port}")
    print('=' * 70)

    start_time = time.time()

    server = EvalServer(test_case, latch_dir)

    async def process_request(path, request_headers):
        if path == "/agent":
            return None
        return (404, [], b"")

    async with websockets.serve(
        server.handle_console_connection,
        "localhost",
        port,
        process_request=process_request
    ):
        print(f"[eval] WebSocket server listening on ws://localhost:{port}/agent")
        print(f"[eval] Connect your console to this URL")

        while not server.is_done():
            await asyncio.sleep(1)

    duration_ms = (time.time() - start_time) * 1000

    notebook_state = {
        "cells": [],
        "cell_status": server.cell_status,
        "cell_last_run_outputs": {},
    }

    test_result = TestResult(
        test_id=test_case.id,
        conversation_history=server.conversation_history,
        notebook_state=notebook_state,
        duration_ms=duration_ms,
    )

    print(f"\n[eval] Eval completed in {duration_ms/1000:.2f}s")
    print(f"[eval] Total conversation turns: {len(server.conversation_history)}")

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

    print(f"\n[eval] Judging result...")
    eval_result = await judge.evaluate(test_case, result)
    result.eval_result = eval_result

    print(f"\n{'=' * 70}")
    print(f"Score: {eval_result.score:.2f}")
    print(f"Passed: {'✓ PASS' if eval_result.passed else '✗ FAIL'}")
    print(f"\nReasoning:")
    print(eval_result.reasoning)
    print(f"\nSuccesses ({len(eval_result.successes)}):")
    for success in eval_result.successes:
        print(f"  ✓ {success}")
    print(f"\nFailures ({len(eval_result.failures)}):")
    for failure in eval_result.failures:
        print(f"  ✗ {failure}")
    print('=' * 70)

    output_file = eval_dir / f"result_{test_case.id}.json"
    with open(output_file, "w") as f:
        json.dump(result.model_dump(), f, indent=2)
    print(f"\n[eval] Result saved to {output_file}")

if __name__ == "__main__":
    asyncio.run(main())
