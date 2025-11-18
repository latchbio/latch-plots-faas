from .base import BinaryGrader, GraderResult

class BooleanCheckGrader(BinaryGrader):
    def evaluate_answer(self, agent_answer: dict, config: dict) -> GraderResult:
        ground_truth = config.get("ground_truth", {})

        metrics = {}
        all_pass = True
        failures = []

        for field, expected_value in ground_truth.items():
            if field not in agent_answer:
                all_pass = False
                failures.append(f"Missing field: {field}")
                continue

            actual_value = agent_answer[field]

            if not isinstance(expected_value, bool):
                all_pass = False
                failures.append(f"{field}: Expected boolean value in ground_truth, got {type(expected_value)}")
                continue

            if not isinstance(actual_value, bool):
                all_pass = False
                failures.append(f"{field}: Expected boolean value from agent, got {type(actual_value)} with value {actual_value}")
                continue

            matches = actual_value == expected_value

            metrics[f"{field}_actual"] = actual_value
            metrics[f"{field}_expected"] = expected_value
            metrics[f"{field}_pass"] = matches

            if not matches:
                all_pass = False
                failures.append(f"{field}: {actual_value} (expected: {expected_value})")

        reasoning = self._format_boolean_reasoning(ground_truth, agent_answer, metrics, failures, all_pass)

        return GraderResult(
            passed=all_pass,
            metrics=metrics,
            reasoning=reasoning,
            agent_answer=agent_answer
        )

    def _format_boolean_reasoning(self, ground_truth, agent_answer, metrics, failures, passed):
        lines = []
        lines.append(f"Boolean Check: {'PASS' if passed else 'FAIL'}")
        lines.append("")

        for field in ground_truth.keys():
            if f"{field}_actual" in metrics:
                actual = metrics[f"{field}_actual"]
                expected = metrics[f"{field}_expected"]
                field_pass = metrics[f"{field}_pass"]
                check = "✓" if field_pass else "✗"
                lines.append(f"- {field}: {actual} (expected: {expected}) {check}")

        if not passed:
            lines.append("")
            lines.append("Failures:")
            for failure in failures:
                lines.append(f"  - {failure}")

        return "\n".join(lines)
