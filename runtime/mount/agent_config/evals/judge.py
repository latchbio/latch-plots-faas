import json
import textwrap

import anthropic

from eval_types import EvalResult, TestCase, TestResult


class LLMJudge:
    def __init__(self, api_key: str, base_url: str, headers: dict):
        self.client = anthropic.AsyncAnthropic(
            api_key=api_key,
            base_url=base_url,
            default_headers=headers
        )

    async def evaluate(self, test_case: TestCase, test_result: TestResult) -> EvalResult:
        conversation_summary = self._format_conversation(test_result.conversation_history)
        notebook_summary = self._format_notebook(test_result.notebook_state)
        globals_summary = self._format_globals(test_result.notebook_state.get("globals", {}))
        reactivity_summary = test_result.notebook_state.get("reactivity", "")

        prompt = textwrap.dedent(f"""
            You are evaluating an AI agent's performance on a task.

            TASK GIVEN TO AGENT:
            {test_case.task}

            DATA PROVIDED:
            {test_case.data_node or "None"}

            EVALUATION CRITERIA:
            {test_case.judge_prompt}

            AGENT'S CONVERSATION HISTORY:
            {conversation_summary}

            FINAL NOTEBOOK STATE:
            {notebook_summary}

            GLOBAL VARIABLES IN KERNEL:
            {globals_summary}

            REACTIVITY/SIGNAL GRAPH:
            {reactivity_summary}

            Evaluate the agent's performance according to the criteria. Be thorough and fair.
            Consider all aspects: the notebook cells, the global variables created, and the reactivity graph.

            Respond ONLY with valid JSON in this exact format (no markdown, no additional text):
            {{
              "score": <number between 0.0 and 1.0>,
              "passed": <boolean>,
              "reasoning": "<detailed reasoning>",
              "successes": ["<success 1>", "<success 2>", ...],
              "failures": ["<failure 1>", "<failure 2>", ...]
            }}
        """).strip()

        response = await self.client.messages.create(
            model="claude-sonnet-4-5-20250929",
            max_tokens=4096,
            messages=[{"role": "user", "content": prompt}]
        )

        response_text = response.content[0].text.strip()
        if response_text.startswith("```json"):
            response_text = response_text[7:]
        if response_text.startswith("```"):
            response_text = response_text[3:]
        if response_text.endswith("```"):
            response_text = response_text[:-3]
        response_text = response_text.strip()

        result_json = json.loads(response_text)
        return EvalResult(**result_json)

    def _format_conversation(self, history: list[dict]) -> str:
        lines = []
        for msg in history:
            role = msg.get("role", "unknown")
            content = msg.get("content", "")

            if isinstance(content, str):
                lines.append(f"{role.upper()}: {content[:500]}")
            elif isinstance(content, list):
                for block in content:
                    if isinstance(block, dict):
                        block_type = block.get("type")
                        if block_type == "text":
                            text = block.get("text", "")[:500]
                            lines.append(f"{role.upper()}: {text}")
                        elif block_type == "tool_use":
                            tool_name = block.get("name", "unknown")
                            lines.append(f"{role.upper()}: [Called tool: {tool_name}]")
                        elif block_type == "tool_result":
                            result = block.get("content", "")[:200]
                            lines.append(f"{role.upper()}: [Tool result: {result}]")

        return "\n".join(lines)

    def _format_notebook(self, notebook_state: dict) -> str:
        context = notebook_state.get("context", {})
        cells = context.get("cells", [])
        lines = [f"Total cells: {len(cells)}\n"]

        for i, cell in enumerate(cells):
            cell_type = cell.get("cell_type", "unknown")
            source = cell.get("source", "")
            status = cell.get("status", "unknown")

            source_preview = source[:200] if isinstance(source, str) else str(source)[:200]
            lines.append(f"\nCell {i} ({cell_type}, status: {status}):")
            lines.append(f"  {source_preview}")

            if cell.get("outputs"):
                lines.append(f"  Outputs: {len(cell['outputs'])} items")

        return "\n".join(lines)

    def _format_globals(self, globals_data: dict) -> str:
        if not globals_data:
            return "No global variables"

        lines = [f"Total global variables: {len(globals_data)}\n"]

        for var_name in sorted(globals_data.keys()):
            var_info = globals_data[var_name]
            if isinstance(var_info, dict):
                var_type = var_info.get("type", "unknown")
                lines.append(f"\n{var_name} ({var_type}):")
                for key, value in var_info.items():
                    if key != "type":
                        value_str = str(value)[:200]
                        lines.append(f"  {key}: {value_str}")
            else:
                var_str = str(var_info)[:200]
                lines.append(f"\n{var_name}: {var_str}")

        return "\n".join(lines)
