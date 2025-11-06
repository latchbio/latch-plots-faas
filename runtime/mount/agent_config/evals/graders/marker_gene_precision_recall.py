from .base import BinaryGrader, GraderResult
from eval_types import TestResult

class MarkerGenePrecisionRecallGrader(BinaryGrader):
    def evaluate(self, test_result: TestResult, config: dict) -> GraderResult:
        canonical_markers = config.get("canonical_markers", [])
        scoring = config.get("scoring", {})
        thresholds = scoring.get("pass_thresholds", {})
        precision_threshold = thresholds.get("precision_at_k", 0.60)
        recall_threshold = thresholds.get("recall_at_k", 0.50)

        agent_answer = self.extract_answer_from_tags(test_result.conversation_history)

        if agent_answer is None:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Agent did not provide answer in <EVAL_ANSWER> tags",
                agent_answer=None
            )

        if "top_marker_genes" not in agent_answer:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Agent answer missing required field: top_marker_genes",
                agent_answer=agent_answer
            )

        predicted_genes = agent_answer["top_marker_genes"]
        k = len(predicted_genes)

        canonical_set = set(gene.lower() for gene in canonical_markers)
        predicted_set = set(gene.lower() for gene in predicted_genes)

        true_positives = canonical_set & predicted_set
        false_positives = predicted_set - canonical_set
        false_negatives = canonical_set - predicted_set

        precision_at_k = len(true_positives) / k if k > 0 else 0.0
        recall_at_k = len(true_positives) / len(canonical_set) if len(canonical_set) > 0 else 0.0

        precision_pass = precision_at_k >= precision_threshold
        recall_pass = recall_at_k >= recall_threshold
        passed = precision_pass and recall_pass

        original_case_map = {gene.lower(): gene for gene in predicted_genes}
        canonical_case_map = {gene.lower(): gene for gene in canonical_markers}

        true_positive_genes = [original_case_map.get(g, canonical_case_map.get(g, g)) for g in true_positives]
        false_positive_genes = [original_case_map.get(g, g) for g in false_positives]
        false_negative_genes = [canonical_case_map.get(g, g) for g in false_negatives]

        metrics = {
            "k": k,
            "precision_at_k": precision_at_k,
            "recall_at_k": recall_at_k,
            "precision_threshold": precision_threshold,
            "recall_threshold": recall_threshold,
            "true_positives": sorted(true_positive_genes),
            "false_positives": sorted(false_positive_genes),
            "false_negatives": sorted(false_negative_genes),
            "num_true_positives": len(true_positives),
            "num_false_positives": len(false_positives),
            "num_false_negatives": len(false_negatives),
            "num_canonical_markers": len(canonical_set),
            "precision_pass": precision_pass,
            "recall_pass": recall_pass,
        }

        reasoning = self._format_precision_recall_reasoning(
            k,
            precision_at_k,
            recall_at_k,
            precision_threshold,
            recall_threshold,
            true_positive_genes,
            false_positive_genes,
            false_negative_genes,
            precision_pass,
            recall_pass,
            passed
        )

        return GraderResult(
            passed=passed,
            metrics=metrics,
            reasoning=reasoning,
            agent_answer=agent_answer
        )

    def _format_precision_recall_reasoning(self, k, precision, recall,
                                          precision_threshold, recall_threshold,
                                          true_positives, false_positives, false_negatives,
                                          precision_pass, recall_pass, passed):
        lines = []
        lines.append(f"Marker Gene Precision/Recall: {'PASS' if passed else 'FAIL'}")
        lines.append("")

        precision_check = "✓" if precision_pass else "✗"
        lines.append(f"Precision@{k}: {precision:.3f} {precision_check} (threshold: {precision_threshold:.3f})")

        recall_check = "✓" if recall_pass else "✗"
        lines.append(f"Recall@{k}: {recall:.3f} {recall_check} (threshold: {recall_threshold:.3f})")

        lines.append("")
        lines.append(f"True Positives ({len(true_positives)}):")
        if true_positives:
            for gene in sorted(true_positives):
                lines.append(f"  ✓ {gene}")
        else:
            lines.append("  None")

        lines.append("")
        lines.append(f"False Positives ({len(false_positives)}):")
        if false_positives:
            for gene in sorted(false_positives):
                lines.append(f"  + {gene}")
        else:
            lines.append("  None")

        lines.append("")
        lines.append(f"False Negatives ({len(false_negatives)}):")
        if false_negatives:
            for gene in sorted(false_negatives):
                lines.append(f"  - {gene}")
        else:
            lines.append("  None")

        lines.append("")
        lines.append(f"Result: {'PASS' if passed else 'FAIL'}")
        if not passed:
            failures = []
            if not precision_pass:
                failures.append(f"Precision {precision:.3f} < {precision_threshold:.3f}")
            if not recall_pass:
                failures.append(f"Recall {recall:.3f} < {recall_threshold:.3f}")
            lines.append(f"Reasons: {'; '.join(failures)}")

        return "\n".join(lines)
