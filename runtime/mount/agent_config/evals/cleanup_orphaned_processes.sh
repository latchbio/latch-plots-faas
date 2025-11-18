#!/bin/bash

echo "========================================"
echo "Orphaned Process Cleanup Script"
echo "========================================"
echo ""

echo "Scanning for orphaned processes..."
echo ""

KERNEL_PIDS=$(ps aux | grep "kernel.py" | grep -v grep | awk '{print $2}')
AGENT_PIDS=$(ps aux | grep "agent.py" | grep -v grep | awk '{print $2}')
EVAL_PIDS=$(ps aux | grep "run_batch_evals.py" | grep -v grep | awk '{print $2}')
EVAL_SERVER_PIDS=$(ps aux | grep "eval_server.py" | grep -v grep | awk '{print $2}')

BATCH_WORKER_PIDS=$(ps -eo pid,pgid,command | grep "run_batch_evals.py" | grep -v grep | awk '{print $2}' | sort -u | while read pgid; do ps -g $pgid -o pid= 2>/dev/null; done | sort -u)

KERNEL_COUNT=$(echo "$KERNEL_PIDS" | grep -c .)
AGENT_COUNT=$(echo "$AGENT_PIDS" | grep -c .)
EVAL_COUNT=$(echo "$EVAL_PIDS" | grep -c .)
EVAL_SERVER_COUNT=$(echo "$EVAL_SERVER_PIDS" | grep -c .)
BATCH_WORKER_COUNT=$(echo "$BATCH_WORKER_PIDS" | grep -c .)

echo "Found:"
echo "  - $KERNEL_COUNT kernel.py processes"
echo "  - $AGENT_COUNT agent.py processes"
echo "  - $EVAL_COUNT run_batch_evals.py processes"
echo "  - $EVAL_SERVER_COUNT eval_server.py processes"
echo "  - $BATCH_WORKER_COUNT processes in run_batch_evals process groups (including workers)"
echo ""

TOTAL=$((KERNEL_COUNT + AGENT_COUNT + EVAL_COUNT + EVAL_SERVER_COUNT + BATCH_WORKER_COUNT))

if [ "$TOTAL" -eq 0 ]; then
    echo "No orphaned processes found. Nothing to clean up."
    exit 0
fi

echo "This will kill all these processes."
read -p "Continue? (y/N) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cleanup cancelled"
    exit 1
fi

echo ""
echo "Killing processes..."

if [ -n "$KERNEL_PIDS" ]; then
    echo "Killing kernel.py processes..."
    echo "$KERNEL_PIDS" | xargs kill -9 2>/dev/null
    echo "  ✓ Killed $KERNEL_COUNT kernel.py processes"
fi

if [ -n "$AGENT_PIDS" ]; then
    echo "Killing agent.py processes..."
    echo "$AGENT_PIDS" | xargs kill -9 2>/dev/null
    echo "  ✓ Killed $AGENT_COUNT agent.py processes"
fi

if [ -n "$EVAL_PIDS" ]; then
    echo "Killing run_batch_evals.py processes..."
    echo "$EVAL_PIDS" | xargs kill -9 2>/dev/null
    echo "  ✓ Killed $EVAL_COUNT run_batch_evals.py processes"
fi

if [ -n "$EVAL_SERVER_PIDS" ]; then
    echo "Killing eval_server.py processes..."
    echo "$EVAL_SERVER_PIDS" | xargs kill -9 2>/dev/null
    echo "  ✓ Killed $EVAL_SERVER_COUNT eval_server.py processes"
fi

if [ -n "$BATCH_WORKER_PIDS" ]; then
    echo "Killing run_batch_evals process groups (main + all workers)..."
    echo "$BATCH_WORKER_PIDS" | xargs kill -9 2>/dev/null
    echo "  ✓ Killed $BATCH_WORKER_COUNT processes from run_batch_evals groups"
fi

echo ""
echo "Cleanup complete!"
echo ""
echo "Remaining eval-related processes:"
ps aux | grep -E "(kernel|agent|eval)" | grep -v grep | grep -v cleanup || echo "  None found"
