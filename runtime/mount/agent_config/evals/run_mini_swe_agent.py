#!/usr/bin/env python3

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from eval_types import TestCase, TestResult
from graders import GRADER_REGISTRY
from minisweagent.agents.default import DefaultAgent
from minisweagent.environments.local import LocalEnvironment
from minisweagent.models import get_model


def run_bioinformatics_eval(eval_file: str, keep_workspace: bool = False):
    print("=" * 80)
    print("Running mini-swe-agent on bioinformatics eval")
    print("=" * 80)

    eval_path = Path(__file__).parent / eval_file
    if not eval_path.exists():
        print(f"Error: Eval file not found: {eval_path}")
        sys.exit(1)

    eval_data = json.loads(eval_path.read_text())
    test_case = TestCase(**eval_data)

    print(f"\nEval: {test_case.id}")
    print("\nTask:")
    print("-" * 80)
    print(test_case.task)
    print("-" * 80)

    cache_dir = Path(tempfile.gettempdir()) / "mini_swe_bio_eval" / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    tmpdir = Path(tempfile.gettempdir()) / "mini_swe_bio_eval"
    work_dir = tmpdir / test_case.id
    if work_dir.exists():
        import shutil
        shutil.rmtree(work_dir)
    work_dir.mkdir()

    print(f"\nWorking directory: {work_dir}")
    if keep_workspace:
        print("   (Workspace will be preserved after run)")
    else:
        print("   (Workspace will be deleted after run)")

    print("\n" + "=" * 80)
    print("Staging data files...")
    print("=" * 80)

    data_nodes = test_case.data_node if isinstance(test_case.data_node, list) else [test_case.data_node]

    original_dir = os.getcwd()

    try:
        os.chdir(str(work_dir))
        print(f"Changed to working directory: {work_dir}")

        contextual_data = []
        for node in data_nodes:
            data_filename = Path(node).name
            cached_file = cache_dir / data_filename

            if cached_file.exists():
                print(f"\nUsing cached data: {data_filename}")
                os.symlink(cached_file, data_filename)
            else:
                print(f"\nDownloading data: {node}")
                subprocess.run(
                    ["latch", "cp", node, str(cached_file)],
                    check=True,
                    capture_output=True
                )
                print(f"Data cached: {cached_file}")
                os.symlink(cached_file, data_filename)

            print(f"Linked: {data_filename} -> {cached_file}")

            contextual_data.append({
                "type": "File",
                "path": node,
                "id": node.replace("latch:///", "").replace(".csv", "").replace(".h5ad", ""),
            })

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

        print("\nInitializing agent...")
        model = get_model()
        env = LocalEnvironment()
        agent = DefaultAgent(model, env)

        print(f"   Model: {model}")
        print(f"   Environment: {env}")

        print("\n" + "=" * 80)
        print("Running agent on task...")
        print("=" * 80)

        import io

        captured_output = io.StringIO()
        original_stdout = sys.stdout
        original_stderr = sys.stderr

        class TeeOutput:
            def __init__(self, *streams):
                self.streams = streams

            def write(self, data):
                for stream in self.streams:
                    stream.write(data)
                    stream.flush()

            def flush(self):
                for stream in self.streams:
                    stream.flush()

        sys.stdout = TeeOutput(original_stdout, captured_output)
        sys.stderr = TeeOutput(original_stderr, captured_output)

        agent_answer = None

        try:
            print("\nSending task to agent...")
            result = agent.run(task_prompt)
            print("\nAgent completed successfully")
            print(f"   Result: {result}")
        except Exception as e:
            if "Submitted" in str(type(e).__name__):
                print("\nAgent finished (raised Submitted exception)")
            else:
                print(f"\nAgent error: {e}")
                import traceback
                traceback.print_exc()
        finally:
            sys.stdout = original_stdout
            sys.stderr = original_stderr

            agent_log_file = work_dir / "agent_output.log"
            agent_log_file.write_text(captured_output.getvalue())
            print(f"\nAgent output saved to: {agent_log_file}")

            trajectory_file = work_dir / "trajectory.json"
            if hasattr(agent, "messages"):
                trajectory_data = {
                    "messages": agent.messages,
                    "actions": getattr(agent, "actions", [])
                }
                trajectory_file.write_text(json.dumps(trajectory_data, indent=2))
                print(f"Agent trajectory saved to: {trajectory_file}")
                print(f"   Total message exchanges: {len(agent.messages)}")

            eval_answer_file = work_dir / "eval_answer.json"
            if eval_answer_file.exists():
                eval_answer_content = eval_answer_file.read_text()
                try:
                    agent_answer = json.loads(eval_answer_content)
                    print(f"Parsed agent answer from eval_answer.json")
                except json.JSONDecodeError as e:
                    print(f"Warning: Failed to parse eval_answer.json: {e}")

    finally:
        os.chdir(original_dir)

    print("\n" + "=" * 80)
    print("Checking results...")
    print("=" * 80)

    print("\nFiles in working directory:")
    for item in sorted(work_dir.iterdir()):
        if item.is_file() and not item.name.startswith("."):
            size = item.stat().st_size
            print(f"   - {item.name} ({size} bytes)")

    trajectory_path = work_dir / "trajectory.json"
    if trajectory_path.exists():
        print("\nAgent trajectory:")
        print("   " + "-" * 76)
        trajectory = json.loads(trajectory_path.read_text())
        messages = trajectory.get("messages", [])

        for i, msg in enumerate(messages):
            role = msg.get("role", "unknown")
            content = msg.get("content", "")

            if i == 0:
                print(f"   [{i}] {role}: {content[:200]}..." if len(content) > 200 else f"   [{i}] {role}: {content}")
            elif "<returncode>" in content:
                returncode_match = re.search(r"<returncode>(\d+)</returncode>", content)
                output_match = re.search(r"<output>(.*?)</output>", content, re.DOTALL)

                returncode = returncode_match.group(1) if returncode_match else "?"
                output = output_match.group(1).strip() if output_match else ""

                print(f"   [{i}] {role}: returncode={returncode}")
                if output:
                    preview = output[:150] + "..." if len(output) > 150 else output
                    print(f"        output: {preview}")
            else:
                bash_match = re.search(r"```bash\n(.*?)\n```", content, re.DOTALL)
                if bash_match:
                    command = bash_match.group(1).strip()
                    print(f"   [{i}] {role}: {command}")
                else:
                    preview = content[:150] + "..." if len(content) > 150 else content
                    print(f"   [{i}] {role}: {preview}")
        print("   " + "-" * 76)

    eval_answer_path = work_dir / "eval_answer.json"
    if eval_answer_path.exists():
        print("\nAgent answer (eval_answer.json):")
        print("   " + "-" * 76)
        content = eval_answer_path.read_text()
        print("   " + content.replace("\n", "\n   "))
        print("   " + "-" * 76)
    else:
        print("\neval_answer.json not found!")

    if test_case.grader and agent_answer is not None:
        print("\n" + "=" * 80)
        print("Running grader...")
        print("=" * 80)

        grader_type = test_case.grader.get("type")
        grader_config = test_case.grader.get("config", {})

        if grader_type in GRADER_REGISTRY:
            grader_cls = GRADER_REGISTRY[grader_type]
            grader = grader_cls()
            grader_result = grader.evaluate_answer(agent_answer, grader_config)

            print(f"\n{'EVAL PASSED' if grader_result.passed else 'EVAL FAILED'}")
            print("\nGrader reasoning:")
            print("-" * 80)
            print(grader_result.reasoning)
            print("-" * 80)

            if grader_result.metrics:
                print("\nMetrics:")
                for key, value in grader_result.metrics.items():
                    print(f"   {key}: {value}")

            if grader_result.agent_answer:
                print("\nAgent answer:")
                print(json.dumps(grader_result.agent_answer, indent=2))
        else:
            print(f"\nWarning: Unknown grader type '{grader_type}'")

    print("\n" + "=" * 80)
    print("Cleanup...")
    print("=" * 80)

    if keep_workspace:
        print(f"\nWorkspace preserved at: {work_dir}")
        print(f"Agent logs: {work_dir}/agent_output.log")
        print(f"Agent trajectory: {work_dir}/trajectory.json")
        print(f"Agent answer: {work_dir}/eval_answer.json")
        print(f"\nTo inspect results:")
        print(f"  cd {work_dir}")
    else:
        import shutil
        shutil.rmtree(work_dir)
        print(f"\nWorkspace deleted: {work_dir}")
        print(f"Data cache preserved at: {cache_dir}")


def main():
    parser = argparse.ArgumentParser(description="Run mini-swe-agent on bioinformatics eval")
    parser.add_argument("--eval", required=True, help="Eval file to run (e.g., curio_seeker_ovary/curio_mt_percentage.json)")
    parser.add_argument("--keep-workspace", action="store_true", help="Keep the workspace directory after completion")
    args = parser.parse_args()

    run_bioinformatics_eval(args.eval, args.keep_workspace)


if __name__ == "__main__":
    main()
