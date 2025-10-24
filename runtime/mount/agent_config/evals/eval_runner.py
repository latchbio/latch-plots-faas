import argparse
import asyncio
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent))

from eval_types import TestCase, TestResult
from judge import LLMJudge


class EvalRunner:
    def __init__(self, test_case: TestCase):
        self.test_case = test_case
        self.conversation_history = []
        self.last_message_time = None
        self.agent_sent_result = False
        self.has_questions = False

    async def run_eval(self) -> TestResult:
        print(f"\n{'=' * 70}")
        print(f"Running eval: {self.test_case.id}")
        print('=' * 70)

        start_time = time.time()

        from run_local_eval import setup_and_run_server

        notebook_state = await setup_and_run_server(
            test_case=self.test_case,
            on_message=self._handle_message,
            is_done=self._is_done,
        )

        duration_ms = (time.time() - start_time) * 1000

        test_result = TestResult(
            test_id=self.test_case.id,
            conversation_history=self.conversation_history,
            notebook_state=notebook_state,
            duration_ms=duration_ms,
        )

        print(f"\n[eval] Eval completed in {duration_ms/1000:.2f}s")
        print(f"[eval] Total conversation turns: {len(self.conversation_history)}")

        return test_result

    async def _handle_message(self, msg: dict):
        self.conversation_history.append(msg)
        self.last_message_time = time.time()

        msg_type = msg.get("type")

        if msg_type == "agent_result":
            self.agent_sent_result = True
            structured = msg.get("structured_output", {})
            questions = structured.get("questions")
            self.has_questions = questions is not None and len(questions) > 0

            print(f"[eval] Agent sent result (questions: {self.has_questions})")

    def _is_done(self, cell_status: dict[str, str]) -> bool:
        if not self.agent_sent_result:
            return False

        if self.has_questions:
            print("[eval] Agent has questions, waiting for user...")
            return False

        running_cells = [cid for cid, status in cell_status.items() if status == "running"]
        if running_cells:
            print(f"[eval] Cells still running: {running_cells}")
            return False

        if self.last_message_time is None:
            return False

        idle_time = time.time() - self.last_message_time
        if idle_time < 10:
            return False

        print(f"[eval] Completion detected: result sent, no questions, no running cells, idle for {idle_time:.1f}s")
        return True


async def run_single_eval(eval_file: Path) -> TestResult:
    with open(eval_file) as f:
        test_data = json.load(f)

    test_case = TestCase(**test_data)
    runner = EvalRunner(test_case)
    result = await runner.run_eval()

    from run_local_eval import get_eval_config
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

    return result


async def run_all_evals(eval_dir: Path) -> list[TestResult]:
    results = []

    for eval_file in sorted(eval_dir.glob("*.json")):
        try:
            result = await run_single_eval(eval_file)
            results.append(result)
        except Exception as e:
            print(f"\n[eval] Error running {eval_file.name}: {e}")
            import traceback
            traceback.print_exc()

    print(f"\n\n{'=' * 70}")
    print("SUMMARY")
    print('=' * 70)

    total = len(results)
    passed = sum(1 for r in results if r.eval_result and r.eval_result.passed)
    failed = total - passed

    avg_score = sum(r.eval_result.score for r in results if r.eval_result) / total if total > 0 else 0

    print(f"Total: {total}")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    print(f"Average Score: {avg_score:.2f}")

    print("\nDetailed Results:")
    for result in results:
        if result.eval_result:
            status = "✓ PASS" if result.eval_result.passed else "✗ FAIL"
            print(f"  {status} {result.test_id} - Score: {result.eval_result.score:.2f}")
        else:
            print(f"  ✗ ERROR {result.test_id}")

    return results


def save_results(results: list[TestResult], output_file: Path):
    output_data = {
        "timestamp": time.time(),
        "results": [r.model_dump() for r in results]
    }

    with open(output_file, "w") as f:
        json.dump(output_data, f, indent=2)

    print(f"\n[eval] Results saved to {output_file}")


async def main():
    parser = argparse.ArgumentParser(description="Run agent evals")
    parser.add_argument("--eval", help="Specific eval file to run")
    parser.add_argument("--all", action="store_true", help="Run all evals")
    parser.add_argument("--output", help="Output file for results", default="results.json")
    args = parser.parse_args()

    eval_dir = Path(__file__).parent

    if args.all:
        results = await run_all_evals(eval_dir)
    elif args.eval:
        eval_file = eval_dir / args.eval
        result = await run_single_eval(eval_file)
        results = [result]
    else:
        print("Please specify --eval <file> or --all")
        return

    output_file = eval_dir / args.output
    save_results(results, output_file)


if __name__ == "__main__":
    asyncio.run(main())
