#!/usr/bin/env python3

import asyncio
import os
import sys
from pathlib import Path


def setup_environment_and_paths():
    sandbox_dir = Path(__file__).parent / "sandbox"
    sandbox_dir.mkdir(exist_ok=True)

    root_dir = sandbox_dir / "root"
    root_dir.mkdir(exist_ok=True)

    latch_dir = root_dir / ".latch"
    latch_dir.mkdir(exist_ok=True)

    user_latch_token = Path.home() / ".latch" / "token"
    if user_latch_token.exists():
        token = user_latch_token.read_text().strip()
    else:
        token = "local-dev-token"

    (latch_dir / "token").write_text(token)
    (latch_dir / "id").write_text("99999")
    (latch_dir / "session-id").write_text("local-session-123")
    (latch_dir / "nucleus-url").write_text("https://nucleus.latch.bio")

    os.environ.update({
        "DD_VERSION": "local-dev",
        "DD_SERVICE": "latch-plots-local",
        "DD_ENV": "local",
        "DD_AGENT_HOST": "localhost",
        "DD_TRACE_ENABLED": "false",
        "DD_PROFILING_ENABLED": "false",
        "DD_RUNTIME_METRICS_ENABLED": "false",
        "OTEL_SDK_DISABLED": "true",
        "auth_jwks_url": "https://example.com/jwks",
        "auth_issuer": "local-dev",
        "auth_audience": "local-dev",
        "auth_self_signed_jwk": "{}",
        "auto_reload": "false",
        "logging_mode": "console",
        "domain": "localhost",
        "AGENT_DEBUG": "1",
        "LATCH_SANDBOX_ROOT": str(latch_dir),
    })

    import pathlib
    original_path_new = pathlib.Path.__new__

    def patched_path_new(cls, *args, **kwargs):
        if args and str(args[0]) == "/root/.latch":
            return original_path_new(cls, str(latch_dir), *args[1:], **kwargs)
        return original_path_new(cls, *args, **kwargs)

    pathlib.Path.__new__ = patched_path_new

    print(f"[Setup] Using sandbox: {sandbox_dir}")
    return sandbox_dir, root_dir, latch_dir


sandbox_dir, root_dir, latch_dir = setup_environment_and_paths()

sys.path.insert(0, str(Path(__file__).parent / "runtime"))
sys.path.insert(0, str(Path(__file__).parent / "runtime/mount"))
sys.path.insert(0, str(Path(__file__).parent / "runtime/mount/python_lib"))

from latch_o11y.o11y import setup as setup_o11y

setup_o11y()


async def mock_gql_query(auth, query, variables=None):
    if "plotsSignerHasNotebookAccess" in query:
        return {"data": {"plotsSignerHasNotebookAccess": True}}
    if "plotsNotebookKernelState" in query:
        return {
            "data": {
                "plotsNotebookKernelState": {
                    "widget_states": [],
                    "cell_output_selections": {},
                    "plot_data_selections": {},
                    "viewer_cell_data": {},
                    "plot_configs": [],
                }
            }
        }
    if "tmpPlotsNotebookKernelSnapshotMode" in query:
        return {"data": {"tmpPlotsNotebookKernelSnapshotMode": False}}
    return {"data": {}}


async def mock_add_pod_event(auth, event_type):
    print(f"[Mock] Pod event: {event_type}")


from runtime.mount import entrypoint, utils

utils.gql_query = mock_gql_query
entrypoint.add_pod_event = mock_add_pod_event


async def run_server():
    from hypercorn.asyncio import serve
    from hypercorn.config import Config as HypercornConfig
    from latch_asgi.server import LatchASGIServer

    from runtime.mount.endpoints import http_routes, websocket_routes
    from runtime.mount.entrypoint import shutdown

    original_start_agent = entrypoint.start_agent_proc

    async def patched_start_agent_proc():
        conn_a = entrypoint.a_proc.conn_a = await entrypoint.SocketIo.from_socket(entrypoint.sock_a)

        print("[run_local] Starting agent subprocess with output capture")
        entrypoint.a_proc.proc = await asyncio.create_subprocess_exec(
            sys.executable,
            "-u",
            (entrypoint.dir_p / "agent.py"),
            str(entrypoint.sock_agent_fd),
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
                print(f"[agent-subprocess] {prefix}{line.decode().rstrip()}", flush=True)

        asyncio.create_task(stream_output(entrypoint.a_proc.proc.stdout, ""))
        asyncio.create_task(stream_output(entrypoint.a_proc.proc.stderr, "[stderr] "))
        print(f"[run_local] Agent subprocess started, PID: {entrypoint.a_proc.proc.pid}", flush=True)

        await conn_a.send({"type": "init"})

        print("[run_local] Waiting for agent ready message...")
        while True:
            msg = await conn_a.recv()
            if msg.get("type") == "ready":
                print("[run_local] Agent ready")
                break

        print("[run_local] Starting message handler")
        entrypoint.async_tasks.append(
            asyncio.create_task(entrypoint.handle_agent_messages(entrypoint.a_proc.conn_a))
        )

    entrypoint.start_agent_proc = patched_start_agent_proc

    latch_server = LatchASGIServer(
        http_routes=http_routes,
        websocket_routes=websocket_routes,
        startup_tasks=[patched_start_agent_proc()],
        shutdown_tasks=[shutdown()],
    )

    cfg = HypercornConfig()
    cfg.bind = ["127.0.0.1:8765"]
    cfg.graceful_timeout = 0.1

    print("\n" + "=" * 70)
    print("  Latch Plots Local Development Server")
    print("=" * 70)
    print("  WebSocket endpoints:")
    print("    - ws://127.0.0.1:8765/run (notebook kernel)")
    print("    - ws://127.0.0.1:8765/agent (AI agent)")
    print()
    print("  HTTP endpoints:")
    print("    - http://127.0.0.1:8765/readyz")
    print()
    print(f"  Sandbox directory: {sandbox_dir}")
    print("  Agent logging: ENABLED (AGENT_DEBUG=1)")
    print()
    print("  Frontend: Set useLocalhost=true to connect")
    print("  Press Ctrl+C to stop")
    print("=" * 70 + "\n")

    shutdown_event = asyncio.Event()

    async def await_shutdown():
        await shutdown_event.wait()

    def shutdown_signal(*args):
        print("\n\nShutting down...")
        shutdown_event.set()

    import signal
    signal.signal(signal.SIGINT, shutdown_signal)
    signal.signal(signal.SIGTERM, shutdown_signal)

    try:
        await serve(latch_server.raw_app, cfg, shutdown_trigger=await_shutdown)
    finally:
        await shutdown()


def main():
    try:
        asyncio.run(run_server())
    except KeyboardInterrupt:
        print("\nStopped")


if __name__ == "__main__":
    main()
