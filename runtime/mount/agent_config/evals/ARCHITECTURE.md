# Eval Harness System - Technical Summary

## System Overview

The eval harness is a system for running automated agent evaluations against a Jupyter-like notebook interface. It allows testing agent behavior by spawning local kernel and agent subprocesses that communicate with a browser-based console frontend.

## Architecture

### Components

1. **Eval Server** (`runtime/mount/agent_config/evals/eval_server.py`)
   - WebSocket server listening on `ws://localhost:8765`
   - Two endpoints: `/run` (kernel) and `/agent` (agent)
   - Spawns kernel and agent as separate subprocesses
   - Routes messages between console frontend ↔ agent subprocess ↔ kernel subprocess

2. **Agent Subprocess** (`runtime/mount/agent.py`)
   - Communicates via Unix socket with eval_server
   - Sends `agent_action` messages for operations (create_cell, delete_cell, run_cell, etc.)
   - Waits for `agent_action_response` messages with 10-second timeout
   - Uses `atomic_operation()` method for all operations

3. **Kernel Subprocess** (`runtime/mount/kernel.py`)
   - Communicates via Unix socket with eval_server
   - Executes Python code in notebook cells
   - Handles special operations: `request_globals_summary`, `request_reactivity_summary`

4. **Console Frontend** (`console/web/src/`)
   - React-based browser interface
   - Connects to eval_server via WebSocket
   - Handles `agent_action` messages via `handleAgentAction()` in `actionHandler.ts`
   - Sends `agent_action_response` messages back

## Message Flow

### Normal Operation Flow
```
Agent → agent_conn socket → eval_server → websocket → Console Frontend
Console Frontend → websocket → eval_server → agent_conn socket → Agent
```

### Kernel-Specific Operations
```
Agent → request_globals_summary/request_reactivity_summary
eval_server intercepts → kernel_conn socket → Kernel
Kernel → kernel_conn socket → eval_server → agent_conn socket → Agent
```

## Key Implementation Details

### 1. EVAL_SANDBOX Mode

When `EVAL_SANDBOX=1` environment variable is set:

**Console WebSocket URLs:**
- Kernel: `ws://localhost:8765/run` (instead of production pod URL)
- Agent: `ws://localhost:8765/agent` (instead of production pod URL)
- Auth token: `"local-dev-token"` (instead of JWT)

**Pod Status Bypass:**
- Console skips all pod status checks (creating, starting, stopping)
- Connection status transitions directly: `disconnected` → `wsConnecting` → `connected`

**Session Owner Override:**
- Run buttons require `isSessionOwner=true`
- In eval mode: `sessionOwner = isEvalSandbox || isSessionOwner`
- Location: `console/web/src/pages/PlotsNotebook/cells/code/index.tsx`

**Files Modified:**
- `console/web/vite.config.ts` - Added EVAL_SANDBOX to define object
- `console/web/src/components/plots/transform/utils.ts` - Pod status bypass
- `console/web/src/components/plots/agent/connection.ts` - Agent WebSocket URL
- `console/web/src/pages/PlotsNotebook/cells/code/index.tsx` - Session owner override

### 2. Agent Operations

**Available Operations:**
- `get_context` - Returns notebook state (cells, outputs, etc.)
- `create_cell` - Creates new code cell
- `create_markdown_cell` - Creates new markdown cell
- `edit_cell` - Edits cell source
- `delete_cell` - Deletes a cell
- `run_cell` - Executes a cell
- `stop_cell` - Stops running cell
- `set_widget` - Updates widget value

**Operation Flow:**
1. Agent calls `await self.atomic_operation(action, params)`
2. Generates tx_id, creates Future, stores in `pending_operations[tx_id]`
3. Sends `{"type": "agent_action", "action": action, "params": params, "tx_id": tx_id}`
4. Waits for response with 10-second timeout
5. Receives `{"type": "agent_action_response", "tx_id": tx_id, "status": "success"|"error", ...}`
6. Future is resolved with response data

### 3. Critical Fix: Async Kernel Request Handling

**Problem:**
The eval_server's forwarding loop was blocking when handling kernel operations:
```python
# OLD (BLOCKING):
if action in ("request_globals_summary", "request_reactivity_summary"):
    await self.kernel_conn.send(msg)
    kernel_response = await self.kernel_conn.recv()  # BLOCKS HERE
    await self.agent_conn.send(kernel_response)
    continue
```

When the agent sent multiple operations concurrently (e.g., `get_context`, `request_globals_summary`, `request_reactivity_summary`), the forwarding loop would:
1. Process `get_context` → forward to console ✓
2. Process `request_globals_summary` → BLOCK waiting for kernel response
3. Any subsequent messages from agent (like `delete_cell`) sit in queue, never processed
4. Agent times out after 10 seconds

**Solution:**
Handle kernel operations asynchronously in background tasks:
```python
# NEW (NON-BLOCKING):
async def handle_kernel_request(action, msg, tx_id):
    await self.kernel_conn.send(msg)
    kernel_response = await self.kernel_conn.recv()
    if kernel_response.get("tx_id") == tx_id:
        await self.agent_conn.send(kernel_response)

# In forwarding loop:
if action in ("request_globals_summary", "request_reactivity_summary"):
    tx_id = msg.get("tx_id")
    asyncio.create_task(handle_kernel_request(action, msg, tx_id))
    continue  # Don't block, continue processing messages
```

This allows the forwarding loop to continue processing new messages while kernel operations complete in the background.

## Data Structures

### Cell Object (from get_context)
```python
{
    "cell_id": "string",  # Required for operations
    "cell_type": "code" | "markdown",
    "source": "string",
    "metadata": {...},
    "language": "python"
}
```

### Agent Action Message
```typescript
{
    type: "agent_action",
    action: "create_cell" | "delete_cell" | "run_cell" | ...,
    params: {...},
    tx_id: "string"
}
```

### Agent Action Response
```typescript
{
    type: "agent_action_response",
    tx_id: "string",
    status: "success" | "error",
    error?: "string",
    // Additional fields depending on action
}
```

## Important Code Locations

### Console Frontend
- **Agent connection**: `console/web/src/components/plots/agent/connection.ts`
- **Agent message handler**: `console/web/src/components/plots/agent/useAgent.ts`
- **Action handler**: `console/web/src/components/plots/agent/actionHandler.ts`
- **Kernel connection state**: `console/web/src/components/plots/transform/utils.ts` (useKernelRuntimeStuff)
- **Code cell container**: `console/web/src/pages/PlotsNotebook/cells/code/index.tsx`
- **Cell container (run buttons)**: `console/web/src/pages/PlotsNotebook/cells/base/editMode/index.tsx`

### Backend
- **Eval server**: `runtime/mount/agent_config/evals/eval_server.py`
- **Agent subprocess**: `runtime/mount/agent.py`
- **Kernel subprocess**: `runtime/mount/kernel.py`
- **Eval runner**: `runtime/mount/agent_config/evals/run_local_eval.py`

## Environment Setup

### Required Environment Variables
- `EVAL_SANDBOX=1` - Enables eval mode in console frontend
- `LATCH_SANDBOX_ROOT` - Path for sandboxed file operations
- `DD_TRACE_ENABLED=false` - Disables DataDog tracing in eval mode

### Vite Configuration
Must add to `console/web/vite.config.ts`:
```typescript
const define = {
  "process.env.EVAL_SANDBOX": JSON.stringify(process.env.EVAL_SANDBOX),
  // ... other env vars
};
```

## Common Pitfalls

1. **Cells without IDs**: If cells don't have `cell_id` field, operations like `delete_cell` will fail silently. Always check that `get_context` returns cells with valid IDs.

2. **Blocking on kernel operations**: Any synchronous wait for kernel responses will block the entire message forwarding loop. Always use `asyncio.create_task()` for kernel operations.

3. **Session owner check**: Run buttons are disabled if `isSessionOwner=false`. Must override this in eval mode.

4. **Agent timeout**: Agent operations timeout after 10 seconds. If eval_server doesn't respond, the operation fails with timeout error.

5. **WebSocket connection timing**: The console must connect to the agent WebSocket and send `init` message before the agent can start processing queries.

## Debugging Tips

1. Check console logs for `[AgentConnection] Connected` to verify WebSocket connection
2. Check eval_server logs for `[eval] agent→console:` to see message flow
3. Check agent logs for `[agent] -> <action>` to see operations being sent
4. Check for timeout errors: `OPERATION FAILED: 'action' timed out after 10 seconds`
5. Verify `EVAL_SANDBOX` is set correctly: `process.env.EVAL_SANDBOX === "1"`

## Test Case Structure

Located in `runtime/mount/agent_config/evals/`:
```python
@dataclass
class TestCase:
    id: str           # Unique test identifier
    task: str         # Task description for agent
    data_node: str | None  # Optional data file path
```

Eval completion is detected by `[EVAL_COMPLETE]` marker in agent's `submit_response` summary.
