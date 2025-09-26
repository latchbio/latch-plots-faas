import asyncio
import os
import socket
import sys
import traceback
from dataclasses import dataclass
from enum import Enum

from lplots import _inject
from socketio_thread import SocketIoThread


class AgentMode(Enum):
    PLANNING = "PLANNING"
    EXECUTING = "EXECUTING"
    DEBUGGING = "DEBUGGING"


@dataclass(kw_only=True)
class Agent:
    conn: SocketIoThread
    initialized: bool = False
    api_key: str | None = None
    mode: AgentMode = AgentMode.PLANNING

    async def send(self, msg: dict[str, object]) -> None:
        await self.conn.send(msg)

    async def handle_init(self, msg: dict[str, object]) -> None:
        print(f"[agent] Initializing")

        if self.api_key is None:
            self.api_key = os.environ.get("OPENAI_API_KEY")

        if self.api_key:
            self.initialized = True
            await self.send({
                "type": "agent_status",
                "status": "ready"
            })
            print("[agent] Initialization complete")
        else:
            await self.send({
                "type": "agent_error",
                "error": "OPENAI_API_KEY not set",
                "fatal": True
            })

    async def handle_query(self, msg: dict[str, object]) -> None:
        request_id = msg.get("request_id")
        query = msg.get("query", "")

        print(f"[agent] Processing query {request_id}: {query[:100]}...")

        try:
            await self.send({
                "type": "agent_result",
                "request_id": request_id,
                "status": "success",
                "response": f"Received query: {query}",
                "mode": self.mode.value
            })
        except Exception as e:
            await self.send({
                "type": "agent_result",
                "request_id": request_id,
                "status": "error",
                "error": str(e)
            })

    async def handle_cancel(self, msg: dict[str, object]) -> None:
        request_id = msg.get("request_id")
        print(f"[agent] Cancel request for {request_id}")

    async def handle_kernel_message(self, msg: dict[str, object]) -> None:
        print(f"[agent] Received kernel message: {msg.get('type')}")

    async def accept(self) -> None:
        msg = await self.conn.recv()
        msg_type = msg.get("type")

        print(f"[agent] Received message type: {msg_type}")

        if msg_type == "init":
            await self.handle_init(msg)
        elif msg_type == "agent_query":
            await self.handle_query(msg)
        elif msg_type == "agent_cancel":
            await self.handle_cancel(msg)
        elif msg_type == "kernel_message":
            await self.handle_kernel_message(msg)
        else:
            print(f"[agent] Unknown message type: {msg_type}")


async def main() -> None:
    global loop
    loop = asyncio.get_running_loop()

    sock = socket.socket(family=socket.AF_UNIX, fileno=int(sys.argv[-1]))
    sock.setblocking(False)

    socket_io_thread = SocketIoThread(socket=sock)
    socket_io_thread.start()
    try:
        socket_io_thread.initialized.wait()

        a = Agent(conn=socket_io_thread)
        _inject.agent = a

        await a.send({"type": "ready"})

        while True:
            try:
                await a.accept()
            except Exception:
                traceback.print_exc()

        print("Kernel shutting down...")
    finally:
        socket_io_thread.shutdown.set()
        socket_io_thread.join()


if __name__ == "__main__":
    if sys.platform == "linux":
        from ctypes import CDLL

        libc = CDLL("libc.so.6")
        PR_SET_NAME = 15  # https://github.com/torvalds/linux/blob/2df0c02dab829dd89360d98a8a1abaa026ef5798/include/uapi/linux/prctl.h#L56
        libc.prctl(PR_SET_NAME, b"agent")

    asyncio.run(main())
