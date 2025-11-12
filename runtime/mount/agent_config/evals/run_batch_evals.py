#!/usr/bin/env python3

import argparse
import asyncio
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import textwrap
import threading
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from eval_types import TestCase
from graders import GRADER_REGISTRY
from minisweagent.agents.default import DefaultAgent
from minisweagent.environments.local import LocalEnvironment
from minisweagent.models import get_model

cache_locks = {}


class ProgressDisplay:
    def __init__(self, total_evals: int):
        self.total = total_evals
        self.passed = 0
        self.failed = 0
        self.current_idx = 0
        self.start_time = time.time()

    def start_eval(self, eval_id: str):
        self.current_idx += 1
        dots = '.' * max(1, 65 - len(eval_id))
        print(f"[{self.current_idx}/{self.total}] {eval_id} {dots} RUNNING", flush=True)

    def complete_eval(self, eval_id: str, passed: bool, duration_s: float, reason: str | None = None):
        status = "✓ PASS" if passed else "✗ FAIL"
        dots = '.' * max(1, 65 - len(eval_id))

        if passed:
            self.passed += 1
        else:
            self.failed += 1

        print(f"[{self.current_idx}/{self.total}] {eval_id} {dots} {status} ({duration_s:.1f}s)")
        if reason:
            print(f"       Reason: {reason[:100]}")

        self.print_progress()

    def print_progress(self):
        elapsed = time.time() - self.start_time
        completed = self.passed + self.failed
        avg_time = elapsed / completed if completed > 0 else 0
        print(f"Progress: {completed}/{self.total} | Pass: {self.passed} | Fail: {self.failed} | Elapsed: {elapsed:.1f}s | Avg: {avg_time:.1f}s/eval\n", flush=True)


class ConcurrentProgressDisplay:
    def __init__(self, total_evals: int, max_workers: int):
        self.total = total_evals
        self.passed = 0
        self.failed = 0
        self.completed = 0
        self.start_time = time.time()
        self.lock = threading.Lock()

    def start_eval(self, eval_id: str) -> int:
        with self.lock:
            timestamp = datetime.now().strftime('%H:%M:%S')
            print(f"[{timestamp}] START {eval_id}")
            return 0

    def update_stage(self, eval_id: str, stage: str):
        with self.lock:
            timestamp = datetime.now().strftime('%H:%M:%S')
            print(f"[{timestamp}] {eval_id} -> {stage}")

    def complete_eval(self, eval_id: str, passed: bool, duration_s: float, reason: str | None = None):
        with self.lock:
            self.completed += 1

            if passed:
                self.passed += 1
            else:
                self.failed += 1

            timestamp = datetime.now().strftime('%H:%M:%S')
            status = "✓ PASS" if passed else "✗ FAIL"

            print(f"[{timestamp}] {status} {eval_id} ({duration_s:.1f}s)")
            if reason:
                print(f"           Reason: {reason[:100]}")

            if self.completed % 5 == 0 or self.completed == self.total:
                self.print_progress()

    def print_progress(self):
        elapsed = time.time() - self.start_time
        avg_time = elapsed / self.completed if self.completed > 0 else 0

        print(f"\nProgress: {self.completed}/{self.total} | Pass: {self.passed} | Fail: {self.failed} | Avg: {avg_time:.1f}s/eval\n")


def discover_evals(eval_dir: Path, pattern: str, filter_id: str | None) -> list[Path]:
    eval_files = []

    for json_file in eval_dir.glob(pattern):
        if not json_file.is_file():
            continue

        try:
            eval_data = json.loads(json_file.read_text())
            eval_id = eval_data.get("id")

            if filter_id and not re.search(filter_id, eval_id):
                continue

            eval_files.append(json_file)
        except Exception as e:
            print(f"Warning: Failed to parse {json_file}: {e}")
            continue

    return sorted(eval_files)


def stage_data_files(test_case: TestCase, work_dir: Path, cache_dir: Path, download_timeout: int) -> list[dict]:
    data_nodes = test_case.data_node if isinstance(test_case.data_node, list) else [test_case.data_node]

    contextual_data = []
    for node in data_nodes:
        data_filename = Path(node).name
        cached_file = cache_dir / data_filename

        if node not in cache_locks:
            cache_locks[node] = threading.Lock()

        with cache_locks[node]:
            if not cached_file.exists():
                try:
                    subprocess.run(
                        ["latch", "cp", node, str(cached_file)],
                        check=True,
                        capture_output=True,
                        timeout=download_timeout
                    )
                except subprocess.TimeoutExpired:
                    raise Exception(f"Download timed out after {download_timeout}s: {node}")

        target_file = work_dir / data_filename
        if target_file.exists():
            target_file.unlink()
        os.symlink(cached_file, target_file)

        contextual_data.append({
            "type": "File",
            "path": node,
            "id": node.replace("latch:///", "").replace(".csv", "").replace(".h5ad", ""),
        })

    return contextual_data


def run_agent_with_timeout(test_case: TestCase, work_dir: Path, contextual_data: list[dict], timeout: int) -> dict | None:
    data_context = f"\n\nHere is the context of the selected nodes the user would like to use: <ContextualNodeData>{json.dumps(contextual_data)}</ContextualNodeData>"

    task_prompt = textwrap.dedent(f"""
        {test_case.task}

        IMPORTANT: When you have completed this task:
        1. Write your final answer as a JSON object to a file named `eval_answer.json`
        2. The file should contain ONLY the JSON object with the required fields
        3. After writing the file, run: `echo COMPLETE_TASK_AND_SUBMIT_FINAL_OUTPUT`

        Example eval_answer.json:
        {{
          "field1": value1,
          "field2": value2
        }}
        {data_context}
    """).strip()

    prompt_file = work_dir / "task_prompt.txt"
    prompt_file.write_text(task_prompt)

    result_file = work_dir / "agent_result.json"

    script = f"""
import sys
import os
import json
import traceback
from pathlib import Path

os.environ['LITELLM_LOCAL_MODEL_COST_MAP'] = 'True'

sys.path.insert(0, str(Path(__file__).parent))

work_dir = Path("{work_dir}")
log_file = work_dir / "agent_conversation.log"

def write_conversation_log(agent):
    with open(log_file, "w") as f:
        for msg in agent.messages:
            role = msg.get("role", "unknown")
            content = msg.get("content", "")
            f.write(f"{'=' * 80}\\n")
            f.write(f"ROLE: {{role.upper()}}\\n")
            f.write(f"{'=' * 80}\\n")
            f.write(content)
            f.write(f"\\n\\n")
        f.flush()

try:
    from minisweagent.agents.default import DefaultAgent
    from minisweagent.environments.local import LocalEnvironment
    from minisweagent.models import get_model

    task_prompt = (work_dir / "task_prompt.txt").read_text()

    model = get_model()
    env = LocalEnvironment()
    agent = DefaultAgent(model, env)

    original_step = agent.step
    def step_with_logging():
        result = original_step()
        write_conversation_log(agent)
        return result
    agent.step = step_with_logging

    agent.run(task_prompt)
    write_conversation_log(agent)

    eval_answer_file = Path("eval_answer.json")
    if eval_answer_file.exists():
        result = json.loads(eval_answer_file.read_text())
        (work_dir / "agent_result.json").write_text(json.dumps({{"answer": result, "status": "success"}}))
    else:
        (work_dir / "agent_result.json").write_text(json.dumps({{"answer": None, "status": "no_answer"}}))
except Exception as e:
    error_msg = f"{{type(e).__name__}}: {{str(e)}}"
    error_traceback = traceback.format_exc()

    if "Submitted" in str(type(e).__name__):
        eval_answer_file = Path("eval_answer.json")
        if eval_answer_file.exists():
            result = json.loads(eval_answer_file.read_text())
            (work_dir / "agent_result.json").write_text(json.dumps({{"answer": result, "status": "success"}}))
        else:
            (work_dir / "agent_result.json").write_text(json.dumps({{"error": error_msg, "traceback": error_traceback, "status": "error"}}))
    else:
        (work_dir / "agent_result.json").write_text(json.dumps({{"error": error_msg, "traceback": error_traceback, "status": "error"}}))
"""

    script_file = work_dir / "run_agent.py"
    script_file.write_text(script)

    try:
        proc_result = subprocess.run(
            [sys.executable, str(script_file)],
            cwd=str(work_dir),
            timeout=timeout,
            capture_output=True,
            text=True
        )

        if result_file.exists():
            result_data = json.loads(result_file.read_text())
            if result_data.get("status") == "success":
                return result_data.get("answer")
            elif result_data.get("status") == "error":
                error_msg = result_data.get("error", "Unknown error")
                if result_data.get("traceback"):
                    error_msg += f"\n\nTraceback:\n{result_data['traceback']}"
                raise Exception(error_msg)

        if proc_result.returncode != 0:
            error_parts = []
            if proc_result.stderr:
                error_parts.append(f"stderr: {proc_result.stderr}")
            if proc_result.stdout:
                error_parts.append(f"stdout: {proc_result.stdout}")
            error_msg = "\n".join(error_parts) if error_parts else "Subprocess failed with no output"
            raise Exception(f"Subprocess exited with code {proc_result.returncode}\n{error_msg}")

        return None

    except subprocess.TimeoutExpired:
        raise Exception(f"Agent execution timed out after {timeout} seconds")
    except Exception as e:
        raise


def run_grader(test_case: TestCase, agent_answer: dict) -> dict | None:
    if not test_case.grader:
        return None

    grader_type = test_case.grader.get("type")
    grader_config = test_case.grader.get("config", {})

    if grader_type not in GRADER_REGISTRY:
        return {
            "passed": False,
            "metrics": {},
            "reasoning": f"Unknown grader type: {grader_type}",
            "agent_answer": agent_answer
        }

    grader_cls = GRADER_REGISTRY[grader_type]
    grader = grader_cls()
    grader_result = grader.evaluate_answer(agent_answer, grader_config)

    return {
        "passed": grader_result.passed,
        "metrics": grader_result.metrics,
        "reasoning": grader_result.reasoning,
        "agent_answer": grader_result.agent_answer
    }


def run_single_eval(eval_file: Path, cache_dir: Path, keep_workspace: bool, download_timeout: int, agent_timeout: int, output_dir: Path, progress: ConcurrentProgressDisplay | None = None) -> dict:
    eval_data = json.loads(eval_file.read_text())
    test_case = TestCase(**eval_data)

    eval_download_timeout = test_case.download_timeout if test_case.download_timeout is not None else (test_case.timeout if test_case.timeout is not None else download_timeout)
    eval_agent_timeout = test_case.agent_timeout if test_case.agent_timeout is not None else (test_case.timeout if test_case.timeout is not None else agent_timeout)

    start_time = time.time()

    work_dir = Path(tempfile.gettempdir()) / "mini_swe_bio_eval" / test_case.id
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True)

    log_file = output_dir / "logs" / f"{test_case.id}.log"

    try:
        if progress:
            progress.update_stage(test_case.id, "downloading")
        contextual_data = stage_data_files(test_case, work_dir, cache_dir, eval_download_timeout)

        if progress:
            progress.update_stage(test_case.id, "running-agent")
        agent_answer = run_agent_with_timeout(test_case, work_dir, contextual_data, eval_agent_timeout)

        conversation_log = work_dir / "agent_conversation.log"
        if conversation_log.exists():
            shutil.copy(conversation_log, log_file)
            conversation_log.unlink()

        grader_result = None
        if test_case.grader and agent_answer:
            if progress:
                progress.update_stage(test_case.id, "grading")
            grader_result = run_grader(test_case, agent_answer)

        duration_s = time.time() - start_time

        passed = grader_result.get("passed", False) if grader_result else True

        if passed and agent_answer:
            (work_dir / "run_agent.py").unlink(missing_ok=True)
            (work_dir / "task_prompt.txt").unlink(missing_ok=True)
            (work_dir / "agent_result.json").unlink(missing_ok=True)

        return {
            "test_id": test_case.id,
            "status": "completed",
            "duration_s": duration_s,
            "grader_result": grader_result,
            "agent_answer": agent_answer,
            "workspace_dir": str(work_dir)
        }

    except Exception as e:
        duration_s = time.time() - start_time

        error_msg = f"Error: {type(e).__name__}: {str(e)}"
        log_file.write_text(error_msg)

        return {
            "test_id": test_case.id,
            "status": "error",
            "duration_s": duration_s,
            "error": str(e),
            "workspace_dir": str(work_dir)
        }

    finally:
        if not keep_workspace and work_dir.exists():
            shutil.rmtree(work_dir)


def save_batch_summary(output_dir: Path, results: list[dict], progress: ProgressDisplay, eval_dir: Path, pattern: str) -> Path:
    run_id = f"batch_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

    durations = [r["duration_s"] for r in results]

    model_name = os.environ.get("MSWEA_MODEL_NAME", "default")

    summary = {
        "run_id": run_id,
        "model": model_name,
        "start_time": datetime.fromtimestamp(progress.start_time).isoformat(),
        "end_time": datetime.now().isoformat(),
        "total_duration_s": time.time() - progress.start_time,
        "eval_dir": str(eval_dir),
        "pattern": pattern,
        "total_evals": len(results),
        "passed": progress.passed,
        "failed": progress.failed,
        "results": results,
        "aggregate_metrics": {
            "avg_duration_s": sum(durations) / len(durations) if durations else 0,
            "median_duration_s": sorted(durations)[len(durations) // 2] if durations else 0,
            "pass_rate": progress.passed / len(results) if results else 0
        }
    }

    summary_file = output_dir / f"batch_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    summary_file.write_text(json.dumps(summary, indent=2))

    return summary_file


def print_final_summary(results: list[dict], output_dir: Path, progress: ProgressDisplay):
    print("\n" + "=" * 80)
    print("BATCH EVAL SUMMARY")
    print("=" * 80)

    total_time = time.time() - progress.start_time

    print(f"Total Evals:     {len(results)}")
    print(f"Passed:          {progress.passed} ({progress.passed/len(results)*100:.1f}%)")
    print(f"Failed:          {progress.failed} ({progress.failed/len(results)*100:.1f}%)")
    print(f"Total Time:      {total_time:.1f}s ({int(total_time//60)}m {int(total_time%60)}s)")
    print(f"Avg Time/Eval:   {total_time/len(results):.1f}s\n")

    passed_results = [r for r in results if r.get("grader_result", {}).get("passed", False)]
    if passed_results:
        print(f"PASSED EVALS ({len(passed_results)}):")
        for r in passed_results[:10]:
            print(f"  ✓ {r['eval_id']} ({r['duration_s']:.1f}s)")
        if len(passed_results) > 10:
            print(f"  ... and {len(passed_results) - 10} more")

    failed_results = [r for r in results if r["status"] == "failed" or (r["status"] == "completed" and not r.get("grader_result", {}).get("passed", False))]
    if failed_results:
        print(f"\nFAILED EVALS ({len(failed_results)}):")
        for r in failed_results:
            print(f"  ✗ {r['eval_id']} ({r['duration_s']:.1f}s)")
            if r.get("failure_reason"):
                print(f"    Reason: {r['failure_reason'][:150]}")
            if r.get("grader_metrics"):
                metrics_str = json.dumps(r['grader_metrics'])[:100]
                print(f"    Metrics: {metrics_str}")

    print(f"\nResults saved to:")
    print(f"  - Summary: {output_dir / 'batch_summary_*.json'}")
    print(f"  - Individual results: {output_dir / 'results/'}")
    print(f"  - Logs: {output_dir / 'logs/'}")

    if failed_results:
        failed_ids = "|".join([r["eval_id"] for r in failed_results])
        print(f"\nTo re-run failed evals:")
        print(f'  python run_batch_evals.py --eval-dir . --filter-id "{failed_ids}"')


async def run_single_eval_async(eval_file: Path, cache_dir: Path, keep_workspace: bool, download_timeout: int, agent_timeout: int, output_dir: Path, semaphore: asyncio.Semaphore, progress: ConcurrentProgressDisplay) -> tuple[str, dict]:
    eval_data = json.loads(eval_file.read_text())
    eval_id = eval_data["id"]

    async with semaphore:
        progress.start_eval(eval_id)

        result = await asyncio.to_thread(
            run_single_eval, eval_file, cache_dir, keep_workspace, download_timeout, agent_timeout, output_dir, progress
        )

        result_file = output_dir / "results" / f"{eval_id}.json"
        await asyncio.to_thread(
            result_file.write_text, json.dumps(result, indent=2)
        )

        passed = result.get("grader_result", {}).get("passed", False) if result["status"] == "completed" else False
        failure_reason = None
        if not passed:
            if result["status"] == "error":
                failure_reason = result.get("error", "Unknown error")
            elif result.get("grader_result"):
                failure_reason = result["grader_result"].get("reasoning", "Grader failed")

        progress.complete_eval(eval_id, passed, result["duration_s"], failure_reason)

        return eval_id, {
            "eval_id": eval_id,
            "status": "passed" if passed else "failed",
            "duration_s": result["duration_s"],
            "grader_passed": passed,
            "grader_metrics": result.get("grader_result", {}).get("metrics", {}),
            "failure_reason": failure_reason,
            "result_file": f"results/{eval_id}.json"
        }


shutdown_requested = False

def signal_handler(signum, frame):
    global shutdown_requested
    if not shutdown_requested:
        shutdown_requested = True
        print("\n\nShutdown requested (Ctrl-C). Waiting for running evals to complete...")
        print("Press Ctrl-C again to force quit.\n")
    else:
        print("\nForce quit!")
        sys.exit(1)

async def run_batch_evals_async(args):
    global shutdown_requested

    signal.signal(signal.SIGINT, signal_handler)

    eval_dir = Path(args.eval_dir)
    if not eval_dir.exists():
        print(f"Error: Eval directory not found: {eval_dir}")
        sys.exit(1)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "results").mkdir(exist_ok=True)
    (output_dir / "logs").mkdir(exist_ok=True)

    eval_files = discover_evals(eval_dir, args.pattern, args.filter_id)

    if not eval_files:
        print(f"No eval files found matching pattern '{args.pattern}' in {eval_dir}")
        sys.exit(1)

    if args.dry_run:
        print(f"Would run {len(eval_files)} evals:")
        for ef in eval_files:
            eval_data = json.loads(ef.read_text())
            print(f"  - {eval_data['id']}")
        return

    model_name = os.environ.get("MSWEA_MODEL_NAME", "default (from config)")
    max_workers = args.max_workers

    print("=" * 80)
    print(f"BATCH EVAL RUN - Started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)
    print(f"Model (env): {model_name}")
    print(f"Found {len(eval_files)} eval files to run")
    print(f"Output directory: {output_dir}")
    print(f"Concurrency: {max_workers} worker(s)")
    try:
        test_model = get_model()
        if hasattr(test_model, 'config') and hasattr(test_model.config, 'model_name'):
            print(f"Model: {test_model.config.model_name}")
        elif hasattr(test_model, 'model_name'):
            print(f"Model: {test_model.model_name}")
        else:
            print(f"Model: {type(test_model).__name__}")
    except Exception as e:
        print(f"⚠️  Warning: Could not verify model: {e}")
    print()

    cache_dir = Path(tempfile.gettempdir()) / "mini_swe_bio_eval" / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    print("Scanning eval files for required data files...")
    unique_data_nodes = set()
    for ef in eval_files:
        eval_data = json.loads(ef.read_text())
        if "data_node" in eval_data and eval_data["data_node"]:
            data_node = eval_data["data_node"]
            if isinstance(data_node, list):
                unique_data_nodes.update(data_node)
            else:
                unique_data_nodes.add(data_node)

    if unique_data_nodes:
        print(f"Found {len(unique_data_nodes)} unique data files to download")
        for i, node in enumerate(unique_data_nodes, 1):
            data_filename = Path(node).name
            cached_file = cache_dir / data_filename

            if cached_file.exists():
                print(f"  [{i}/{len(unique_data_nodes)}] {data_filename} - already cached")
            else:
                print(f"  [{i}/{len(unique_data_nodes)}] {data_filename} - downloading...")
                try:
                    subprocess.run(
                        ["latch", "cp", node, str(cached_file)],
                        check=True,
                        capture_output=True
                    )
                    file_size = cached_file.stat().st_size / (1024**3)
                    print(f"      Downloaded {file_size:.2f} GB")
                except Exception as e:
                    print(f"      ⚠️  Download failed: {e}")
        print()
    else:
        print("No data files to pre-download\n")

    if max_workers == 1:
        progress = ProgressDisplay(len(eval_files))
    else:
        progress = ConcurrentProgressDisplay(len(eval_files), max_workers)

    semaphore = asyncio.Semaphore(max_workers)

    if max_workers == 1:
        batch_results = []
        for eval_file in eval_files:
            if shutdown_requested:
                print(f"\nStopping after {len(batch_results)} evals completed")
                break

            eval_data = json.loads(eval_file.read_text())
            eval_id = eval_data["id"]

            progress.start_eval(eval_id)

            result = await asyncio.to_thread(
                run_single_eval, eval_file, cache_dir, args.keep_workspace, args.download_timeout, args.agent_timeout, output_dir
            )

            result_file = output_dir / "results" / f"{eval_id}.json"
            result_file.write_text(json.dumps(result, indent=2))

            passed = result.get("grader_result", {}).get("passed", False) if result["status"] == "completed" else False
            failure_reason = None
            if not passed:
                if result["status"] == "error":
                    failure_reason = result.get("error", "Unknown error")
                elif result.get("grader_result"):
                    failure_reason = result["grader_result"].get("reasoning", "Grader failed")[:100]

            progress.complete_eval(eval_id, passed, result["duration_s"], failure_reason)

            batch_results.append({
                "eval_id": eval_id,
                "status": "passed" if passed else "failed",
                "duration_s": result["duration_s"],
                "grader_passed": passed,
                "grader_metrics": result.get("grader_result", {}).get("metrics", {}),
                "failure_reason": failure_reason,
                "result_file": f"results/{eval_id}.json"
            })
    else:
        tasks = [
            run_single_eval_async(ef, cache_dir, args.keep_workspace, args.download_timeout, args.agent_timeout, output_dir, semaphore, progress)
            for ef in eval_files
        ]

        try:
            results = await asyncio.gather(*tasks, return_exceptions=True)
        except asyncio.CancelledError:
            print("\nTasks cancelled by shutdown request")
            results = []

        batch_results = []
        for result in results:
            if isinstance(result, Exception):
                if not isinstance(result, asyncio.CancelledError):
                    print(f"Error: {result}")
                continue
            _, result_dict = result
            batch_results.append(result_dict)

        if shutdown_requested:
            print(f"\nShutdown complete. Saved results for {len(batch_results)} completed evals.")

    if batch_results:
        save_batch_summary(output_dir, batch_results, progress, eval_dir, args.pattern)
        print_final_summary(batch_results, output_dir, progress)
    else:
        print("\nNo results to save (all evals were interrupted or failed)")


def main():
    parser = argparse.ArgumentParser(description="Run batch evals using mini-swe-agent")
    parser.add_argument("--eval-dir", required=True, help="Directory containing eval JSON files")
    parser.add_argument("--pattern", default="**/*.json", help="Glob pattern to match eval files (default: **/*.json)")
    parser.add_argument("--keep-workspace", action="store_true", help="Keep workspace directories after completion")
    parser.add_argument("--output-dir", default="./batch_results", help="Results output directory (default: ./batch_results)")
    parser.add_argument("--timeout", type=int, default=1200, help="DEPRECATED: Use --agent-timeout instead (default: 1200)")
    parser.add_argument("--download-timeout", type=int, default=600, help="Download timeout in seconds (default: 600)")
    parser.add_argument("--agent-timeout", type=int, default=1200, help="Agent execution timeout in seconds (default: 1200)")
    parser.add_argument("--filter-id", help="Regex pattern to filter eval IDs")
    parser.add_argument("--dry-run", action="store_true", help="Show which evals would run without executing")
    parser.add_argument("--max-workers", type=int, default=1, help="Maximum concurrent evals (default: 1 for sequential)")
    args = parser.parse_args()

    asyncio.run(run_batch_evals_async(args))


if __name__ == "__main__":
    main()
