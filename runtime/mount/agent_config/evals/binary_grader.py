import json
import re
import numpy as np
from dataclasses import dataclass
from sklearn.metrics import roc_auc_score

from eval_types import TestResult

@dataclass
class GraderResult:
    passed: bool
    metrics: dict
    reasoning: str
    agent_answer: dict | None

class BinaryGrader:
    def evaluate(self, test_result: TestResult, config: dict) -> GraderResult:
        raise NotImplementedError

    def extract_answer_from_tags(self, conversation: list[dict]) -> dict | None:
        for msg in reversed(conversation):
            if msg.get("type") != "anthropic_message" or msg.get("role") != "assistant":
                continue

            content = msg.get("content", [])
            for block in content:
                if isinstance(block, dict) and block.get("type") == "tool_use":
                    if block.get("name") == "submit_response":
                        tool_input = block.get("input", {})
                        summary = tool_input.get("summary", "")

                        match = re.search(r'<EVAL_ANSWER>(.*?)</EVAL_ANSWER>', summary, re.DOTALL)
                        if match:
                            json_str = match.group(1).strip()
                            try:
                                return json.loads(json_str)
                            except json.JSONDecodeError as e:
                                print(f"[grader] Failed to parse JSON from EVAL_ANSWER tags: {e}")
                                return None
        return None

class CellTypingGrader(BinaryGrader):
    def evaluate(self, test_result: TestResult, config: dict) -> GraderResult:
        markers = config.get("markers", [])
        thresholds = config.get("thresholds", {
            "auroc": 0.75,
            "gap_ratio": 1.15
        })

        agent_answer = self.extract_answer_from_tags(test_result.conversation_history)

        if agent_answer is None:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Agent did not provide answer in <EVAL_ANSWER> tags",
                agent_answer=None
            )

        required_fields = ["chosen_cluster", "samples_chosen", "samples_other", "cluster_medians"]
        missing = [f for f in required_fields if f not in agent_answer]
        if missing:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning=f"Agent answer missing required fields: {missing}",
                agent_answer=agent_answer
            )

        chosen_cluster = agent_answer["chosen_cluster"]
        samples_chosen = agent_answer["samples_chosen"]
        samples_other = agent_answer["samples_other"]
        cluster_medians = agent_answer["cluster_medians"]

        geneset_scores_chosen = [self._compute_geneset_score(s, markers) for s in samples_chosen]
        geneset_scores_other = [self._compute_geneset_score(s, markers) for s in samples_other]

        try:
            labels = [1] * len(geneset_scores_chosen) + [0] * len(geneset_scores_other)
            scores = geneset_scores_chosen + geneset_scores_other
            auroc = roc_auc_score(labels, scores)
        except Exception as e:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning=f"Failed to compute AUROC: {e}",
                agent_answer=agent_answer
            )

        cluster_geneset_medians = {}
        for cluster_id, marker_vals in cluster_medians.items():
            geneset_score = self._compute_geneset_score(marker_vals, markers)
            cluster_geneset_medians[cluster_id] = geneset_score

        sorted_clusters = sorted(cluster_geneset_medians.items(), key=lambda x: x[1], reverse=True)

        if len(sorted_clusters) < 2:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Need at least 2 clusters to compute gap ratio",
                agent_answer=agent_answer
            )

        top_cluster, top_score = sorted_clusters[0]
        second_cluster, second_score = sorted_clusters[1]

        is_top_cluster = (chosen_cluster == top_cluster)
        chosen_median_score = cluster_geneset_medians.get(chosen_cluster, 0.0)

        if second_score > 0:
            gap_ratio = top_score / second_score
        else:
            gap_ratio = float('inf')

        auroc_pass = auroc >= thresholds["auroc"]
        gap_pass = gap_ratio >= thresholds["gap_ratio"]

        passed = auroc_pass and gap_pass and is_top_cluster

        metrics = {
            "geneset_auroc": auroc,
            "chosen_cluster_median_score": chosen_median_score,
            "top_cluster": top_cluster,
            "top_cluster_median_score": top_score,
            "second_best_cluster": second_cluster,
            "second_best_median_score": second_score,
            "gap_ratio": gap_ratio,
            "is_top_cluster": is_top_cluster,
            "auroc_pass": auroc_pass,
            "gap_pass": gap_pass,
        }

        reasoning = self._format_reasoning(
            chosen_cluster,
            auroc,
            auroc_pass,
            top_cluster,
            top_score,
            second_cluster,
            second_score,
            gap_ratio,
            gap_pass,
            is_top_cluster,
            thresholds,
            passed
        )

        return GraderResult(
            passed=passed,
            metrics=metrics,
            reasoning=reasoning,
            agent_answer=agent_answer
        )

    def _compute_geneset_score(self, sample: dict, markers: list[str]) -> float:
        values = [sample.get(marker, 0.0) for marker in markers]
        if not values:
            return 0.0
        return np.mean(values)

    def _format_reasoning(self, chosen_cluster, auroc, auroc_pass, top_cluster, top_score,
                         second_cluster, second_score, gap_ratio, gap_pass, is_top_cluster,
                         thresholds, passed):
        lines = []
        lines.append(f"Chosen cluster: '{chosen_cluster}' - {'PASS' if passed else 'FAIL'}")
        lines.append("")

        lines.append("Gene-Set Metrics:")
        auroc_check = "✓" if auroc_pass else "✗"
        lines.append(f"- AUROC (one-vs-rest): {auroc:.3f} {auroc_check} (threshold: {thresholds['auroc']})")

        lines.append("")
        lines.append("Cluster Rankings (by median gene-set score):")
        lines.append(f"  1. {top_cluster}: {top_score:.4f}")
        lines.append(f"  2. {second_cluster}: {second_score:.4f}")

        gap_check = "✓" if gap_pass else "✗"
        lines.append(f"  Gap ratio (top/second): {gap_ratio:.3f} {gap_check} (threshold: {thresholds['gap_ratio']})")

        lines.append("")
        if not is_top_cluster:
            lines.append(f"✗ Chosen cluster '{chosen_cluster}' is NOT the top-ranked cluster ('{top_cluster}' is)")
        else:
            lines.append(f"✓ Chosen cluster '{chosen_cluster}' is the top-ranked cluster")

        lines.append("")
        lines.append(f"Result: {'PASS' if passed else 'FAIL'}")
        if not passed:
            failures = []
            if not auroc_pass:
                failures.append("AUROC below threshold")
            if not gap_pass:
                failures.append("Gap ratio below threshold")
            if not is_top_cluster:
                failures.append("Not the top cluster")
            lines.append(f"Reasons: {', '.join(failures)}")

        return "\n".join(lines)

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

class LabelSetJaccardGrader(BinaryGrader):
    def evaluate(self, test_result: TestResult, config: dict) -> GraderResult:
        ground_truth_labels = set(config.get("ground_truth_labels", []))
        scoring = config.get("scoring", {})
        pass_threshold = scoring.get("pass_threshold", 0.90)

        agent_answer = self.extract_answer_from_tags(test_result.conversation_history)

        if agent_answer is None:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Agent did not provide answer in <EVAL_ANSWER> tags",
                agent_answer=None
            )

        if "cell_types_predicted" not in agent_answer:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Agent answer missing required field: cell_types_predicted",
                agent_answer=agent_answer
            )

        predicted_labels = set(agent_answer["cell_types_predicted"])

        intersection = ground_truth_labels & predicted_labels
        union = ground_truth_labels | predicted_labels

        jaccard_index = len(intersection) / len(union) if len(union) > 0 else 0.0

        passed = jaccard_index >= pass_threshold

        true_positives = intersection
        false_positives = predicted_labels - ground_truth_labels
        false_negatives = ground_truth_labels - predicted_labels

        metrics = {
            "jaccard_index": jaccard_index,
            "pass_threshold": pass_threshold,
            "true_positives": sorted(list(true_positives)),
            "false_positives": sorted(list(false_positives)),
            "false_negatives": sorted(list(false_negatives)),
            "predicted_count": len(predicted_labels),
            "ground_truth_count": len(ground_truth_labels),
        }

        reasoning = self._format_jaccard_reasoning(
            jaccard_index,
            pass_threshold,
            true_positives,
            false_positives,
            false_negatives,
            passed
        )

        return GraderResult(
            passed=passed,
            metrics=metrics,
            reasoning=reasoning,
            agent_answer=agent_answer
        )

    def _format_jaccard_reasoning(self, jaccard_index, threshold, true_positives, false_positives, false_negatives, passed):
        lines = []
        lines.append(f"Label Set Comparison: {'PASS' if passed else 'FAIL'}")
        lines.append("")
        lines.append(f"Jaccard Index: {jaccard_index:.3f} (threshold: {threshold:.3f}) {'✓' if jaccard_index >= threshold else '✗'}")
        lines.append("")

        if true_positives:
            lines.append(f"Correct Labels ({len(true_positives)}):")
            for label in sorted(true_positives):
                lines.append(f"  ✓ {label}")
        else:
            lines.append("Correct Labels: None")

        lines.append("")

        if false_positives:
            lines.append(f"Extra Labels ({len(false_positives)}):")
            for label in sorted(false_positives):
                lines.append(f"  + {label}")
        else:
            lines.append("Extra Labels: None")

        lines.append("")

        if false_negatives:
            lines.append(f"Missing Labels ({len(false_negatives)}):")
            for label in sorted(false_negatives):
                lines.append(f"  - {label}")
        else:
            lines.append("Missing Labels: None")

        lines.append("")
        lines.append(f"Result: {'PASS' if passed else 'FAIL'}")

        return "\n".join(lines)

GRADER_REGISTRY = {
    "cell_typing": CellTypingGrader,
    "numeric_tolerance": NumericToleranceGrader,
    "label_set_jaccard": LabelSetJaccardGrader,
}
