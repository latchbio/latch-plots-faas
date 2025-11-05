from .base import BinaryGrader, GraderResult
from eval_types import TestResult

class NumericToleranceGrader(BinaryGrader):
    def evaluate(self, test_result: TestResult, config: dict) -> GraderResult:
        ground_truth = config.get("ground_truth", {})
        tolerances = config.get("tolerances", {})

        agent_answer = self.extract_answer_from_tags(test_result.conversation_history)

        if agent_answer is None:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Agent did not provide answer in <EVAL_ANSWER> tags",
                agent_answer=None
            )

        metrics = {}
        all_pass = True
        failures = []

        for field, expected_value in ground_truth.items():
            if field not in agent_answer:
                all_pass = False
                failures.append(f"Missing field: {field}")
                continue

            actual_value = agent_answer[field]
            tolerance_config = tolerances.get(field, {"type": "absolute", "value": 0})
            tolerance_type = tolerance_config.get("type", "absolute")
            tolerance_value = tolerance_config.get("value", 0)

            if tolerance_type == "absolute":
                within_tolerance = abs(actual_value - expected_value) <= tolerance_value
                error = abs(actual_value - expected_value)
            elif tolerance_type == "relative":
                relative_error = abs(actual_value - expected_value) / abs(expected_value) if expected_value != 0 else float('inf')
                within_tolerance = relative_error <= tolerance_value
                error = relative_error
            else:
                within_tolerance = False
                error = float('inf')

            metrics[f"{field}_actual"] = actual_value
            metrics[f"{field}_expected"] = expected_value
            metrics[f"{field}_error"] = error
            metrics[f"{field}_pass"] = within_tolerance

            if not within_tolerance:
                all_pass = False
                failures.append(f"{field}: {actual_value} vs {expected_value} (error: {error:.2f}, tolerance: {tolerance_value})")

        reasoning = self._format_numeric_reasoning(ground_truth, agent_answer, metrics, failures, all_pass)

        return GraderResult(
            passed=all_pass,
            metrics=metrics,
            reasoning=reasoning,
            agent_answer=agent_answer
        )

    def _format_numeric_reasoning(self, ground_truth, agent_answer, metrics, failures, passed):
        lines = []
        lines.append(f"Numeric Tolerance Check: {'PASS' if passed else 'FAIL'}")
        lines.append("")

        for field in ground_truth.keys():
            if f"{field}_actual" in metrics:
                actual = metrics[f"{field}_actual"]
                expected = metrics[f"{field}_expected"]
                error = metrics[f"{field}_error"]
                field_pass = metrics[f"{field}_pass"]
                check = "✓" if field_pass else "✗"
                lines.append(f"- {field}: {actual} vs {expected} (error: {error:.2f}) {check}")

        if not passed:
            lines.append("")
            lines.append("Failures:")
            for failure in failures:
                lines.append(f"  - {failure}")

        return "\n".join(lines)
