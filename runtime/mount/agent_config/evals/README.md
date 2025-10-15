# Agent Evals

This directory contains end-to-end evaluation tasks for the notebook agent.

## Test Case Format

Each eval is a single JSON file with the following structure:

```json
{
  "id": "unique_test_id",
  "task": "The task to give to the agent",
  "data_node": "latch:///path/to/data.csv (optional)",
  "judge_prompt": "Instructions for the LLM judge on how to evaluate the result"
}
```

## Running Evals

### Step 1: Start the eval server

```bash
cd runtime/mount/agent_config/evals
python eval_server.py --eval cluster_workflow.json
```

The eval server:
- Listens on port 8765 for WebSocket connections (default, use --port to change)
- Spawns its own agent process with the test task
- Acts as a bridge between agent and console

### Step 2: Connect console to eval server

Point your browser console to `ws://localhost:8765/agent`.

The console will communicate with its own kernel (wherever that's running), and the agent will interact with it through the eval server bridge.

## How It Works

1. **Eval Server**: Spawns agent process and bridges to console via WebSocket
2. **Console Connection**: Console connects to eval server and manages its own kernel connection
3. **Auto-Start**: Agent receives initial query via `--initial-query` flag
4. **Real Execution**: Agent's tool calls → eval server → console → kernel
5. **Monitoring**: Eval tracks conversation and completion
6. **Completion Detection**: Eval completes when:
   - Agent sent final result
   - No questions pending
   - No cells running
   - Idle for 10 seconds
7. **LLM Judge**: Claude evaluates the result with structured output
8. **Results**: Score (0-1), pass/fail, successes, failures, reasoning

## Results Format

Results are saved as JSON with:
- Full conversation history
- Cell execution status
- Judge evaluation (score, passed, reasoning, successes, failures)
- Duration and metadata
