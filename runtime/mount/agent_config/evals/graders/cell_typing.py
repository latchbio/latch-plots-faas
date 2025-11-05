import numpy as np
from sklearn.metrics import roc_auc_score

from .base import BinaryGrader, GraderResult
from eval_types import TestResult

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
