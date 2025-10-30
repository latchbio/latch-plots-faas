import asyncio
import os
import sys
import tempfile
from collections.abc import Callable
from pathlib import Path

from eval_types import TestCase


def get_eval_config() -> tuple[str, str, int]:
    user_latch_token = Path.home() / ".latch" / "token"
    if user_latch_token.exists():
        token = user_latch_token.read_text().strip()
    else:
        token = "local-dev-token"

    auth_token_sdk = f"Latch-SDK-Token {token}"
    nucleus_url = "https://nucleus.latch.bio"
    pod_id = 99999

    return auth_token_sdk, nucleus_url, pod_id

async def setup_and_run_server(
    test_case: TestCase,
    on_message: Callable,
    is_done: Callable[[dict[str, str]], bool],
) -> dict:
    with tempfile.TemporaryDirectory() as tmpdir:
        sandbox_dir = Path(tmpdir) / "sandbox"
        sandbox_dir.mkdir()

        latch_dir = sandbox_dir / ".latch"
        latch_dir.mkdir()

        user_latch_token = Path.home() / ".latch" / "token"
        if user_latch_token.exists():
            token = user_latch_token.read_text().strip()
        else:
            token = "local-dev-token"

        (latch_dir / "token").write_text(token)
        (latch_dir / "id").write_text("99999")
        (latch_dir / "session-id").write_text(f"eval-session-{test_case.id}")
        (latch_dir / "nucleus-url").write_text("https://nucleus.latch.bio")

        project_root = Path(__file__).parent.parent.parent.parent.parent
        sys.path.insert(0, str(project_root))
        sys.path.insert(0, str(project_root / "runtime"))
        sys.path.insert(0, str(project_root / "runtime" / "mount"))
        sys.path.insert(0, str(project_root / "runtime" / "mount" / "python_lib"))

        os.environ.update({
            "DD_VERSION": "eval",
            "DD_SERVICE": "latch-plots-eval",
            "DD_ENV": "eval",
            "DD_AGENT_HOST": "localhost",
            "DD_TRACE_ENABLED": "false",
            "DD_PROFILING_ENABLED": "false",
            "DD_RUNTIME_METRICS_ENABLED": "false",
            "OTEL_SDK_DISABLED": "true",
            "auth_jwks_url": "https://example.com/jwks",
            "auth_issuer": "eval",
            "auth_audience": "eval",
            "auth_self_signed_jwk": "{}",
            "auto_reload": "false",
            "logging_mode": "console",
            "domain": "localhost",
            "AGENT_DEBUG": "1",
            "LATCH_SANDBOX_ROOT": str(latch_dir),
        })

        from latch_o11y.o11y import setup as setup_o11y
        setup_o11y()

        from mount import entrypoint, utils

        entrypoint.latch_p = latch_dir
        entrypoint.sdk_token = token
        entrypoint.auth_token_sdk = f"Latch-SDK-Token {token}"
        entrypoint.pod_id = 99999
        entrypoint.pod_session_id = f"eval-session-{test_case.id}"
        utils.auth_token_sdk = f"Latch-SDK-Token {token}"
        utils.nucleus_url = "https://nucleus.latch.bio"
        utils.pod_id = 99999

        async def mock_gql_query(auth, query, variables=None):
            if "plotsSignerHasNotebookAccess" in query:
                return {"data": {"plotsSignerHasNotebookAccess": True}}
            if "plotsNotebookKernelState" in query:
                return {
                    "data": {
                        "plotsNotebookKernelState": {
                            "widget_states": {},
                            "cell_output_selections": {},
                            "plot_data_selections": {},
                            "viewer_cell_data": {},
                            "plot_configs": {},
                        }
                    }
                }
            if "tmpPlotsNotebookKernelSnapshotMode" in query:
                return {"data": {"tmpPlotsNotebookKernelSnapshotMode": False}}
            if "mutation" in query.lower():
                return {"data": {"clientMutationId": "mock"}}
            return {"data": {}}

        async def mock_add_pod_event(auth, event_type):
            pass

        utils.gql_query = mock_gql_query
        entrypoint.gql_query = mock_gql_query
        entrypoint.add_pod_event = mock_add_pod_event

        from runtime.mount.entrypoint import shutdown

        async def patched_start_agent_proc():
            conn_a = entrypoint.a_proc.conn_a = await entrypoint.SocketIo.from_socket(entrypoint.sock_a)

            entrypoint.a_proc.proc = await asyncio.create_subprocess_exec(
                sys.executable,
                "-u",
                (entrypoint.dir_p / "agent.py"),
                str(entrypoint.sock_agent_fd),
                "--initial-query",
                test_case.task,
                pass_fds=[entrypoint.sock_agent_fd],
                stdin=asyncio.subprocess.DEVNULL,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                env={**os.environ, "LATCH_SANDBOX_ROOT": str(latch_dir), "PYTHONUNBUFFERED": "1"},
                preexec_fn=lambda: os.nice(5),
            )

            async def stream_output(stream, prefix=""):
                while True:
                    line = await stream.readline()
                    if not line:
                        break
                    print(f"[agent] {prefix}{line.decode().rstrip()}", flush=True)

            asyncio.create_task(stream_output(entrypoint.a_proc.proc.stdout, ""))
            asyncio.create_task(stream_output(entrypoint.a_proc.proc.stderr, "[stderr] "))

            await conn_a.send({"type": "init"})

            while True:
                msg = await conn_a.recv()
                if msg.get("type") == "ready":
                    break

            async def handle_agent_messages_with_callback():
                while True:
                    msg = await conn_a.recv()
                    await on_message(msg)

                    if entrypoint.current_agent_ctx is None:
                        continue

                    try:
                        import orjson
                        await entrypoint.current_agent_ctx.send_message(orjson.dumps(msg).decode())
                    except Exception as e:
                        print(f"[entrypoint] Error forwarding message: {e}")

            entrypoint.async_tasks.append(
                asyncio.create_task(handle_agent_messages_with_callback())
            )

        entrypoint.start_agent_proc = patched_start_agent_proc

        await entrypoint.start_kernel_proc()
        await patched_start_agent_proc()

        await entrypoint.ready_ev.wait()

        while not is_done(entrypoint.cell_status):
            await asyncio.sleep(1)

        notebook_state = {
            "cells": [],
            "cell_status": dict(entrypoint.cell_status),
            "cell_last_run_outputs": dict(entrypoint.cell_last_run_outputs),
        }

        await shutdown()

        return notebook_state
