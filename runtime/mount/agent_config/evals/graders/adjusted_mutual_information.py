from typing import List

import numpy as np
from sklearn.metrics import adjusted_mutual_info_score

from .base import BinaryGrader, GraderResult
from eval_types import TestResult


class AdjustedMutualInformationGrader(BinaryGrader):
    """Evaluate cluster-condition association via adjusted mutual information."""

    def evaluate(self, test_result: TestResult, config: dict) -> GraderResult:
        scoring = config.get("scoring", {})
        thresholds = scoring.get("pass_thresholds", {})
        ami_min = thresholds.get("ami_min", 0.0)
        ami_max = thresholds.get("ami_max")  # optional upper bound

        agent_answer = self.extract_answer_from_tags(test_result.conversation_history)

        if agent_answer is None:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Agent did not provide answer in <EVAL_ANSWER> tags",
                agent_answer=None,
            )

        required_fields = ["clusters", "conditions", "contingency"]
        missing = [field for field in required_fields if field not in agent_answer]
        if missing:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning=f"Agent answer missing required fields: {missing}",
                agent_answer=agent_answer,
            )

        clusters = agent_answer["clusters"]
        conditions = agent_answer["conditions"]
        contingency = agent_answer["contingency"]

        if not isinstance(clusters, list) or not clusters:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="'clusters' must be a non-empty list",
                agent_answer=agent_answer,
            )

        if not isinstance(conditions, list) or not conditions:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="'conditions' must be a non-empty list",
                agent_answer=agent_answer,
            )

        if not isinstance(contingency, list) or not contingency:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="'contingency' must be a non-empty 2D list",
                agent_answer=agent_answer,
            )

        try:
            matrix = np.array(contingency, dtype=float)
        except Exception as exc:  # pragma: no cover - defensive
            return GraderResult(
                passed=False,
                metrics={},
                reasoning=f"Failed to parse contingency matrix: {exc}",
                agent_answer=agent_answer,
            )

        expected_shape = (len(clusters), len(conditions))
        if matrix.shape != expected_shape:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning=(
                    "Contingency matrix shape does not match provided labels: "
                    f"expected {expected_shape}, got {matrix.shape}"
                ),
                agent_answer=agent_answer,
            )

        if (matrix < 0).any():
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Contingency matrix contains negative counts",
                agent_answer=agent_answer,
            )

        if matrix.sum() <= 0:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Contingency matrix must contain positive counts",
                agent_answer=agent_answer,
            )

        cluster_labels: List[str] = []
        condition_labels: List[str] = []
        for ci, cluster in enumerate(clusters):
            for cj, condition in enumerate(conditions):
                count = matrix[ci, cj]
                if count <= 0:
                    continue
                n = int(round(count))
                if n <= 0:
                    continue
                cluster_labels.extend([cluster] * n)
                condition_labels.extend([condition] * n)

        if not cluster_labels or not condition_labels:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Expanded contingency produced empty label arrays",
                agent_answer=agent_answer,
            )

        ami_value = adjusted_mutual_info_score(cluster_labels, condition_labels)

        passes_min = ami_value >= ami_min
        passes_max = True if ami_max is None else ami_value <= ami_max
        passed = passes_min and passes_max

        metrics = {
            "ami": float(ami_value),
            "ami_min_threshold": ami_min,
            "ami_max_threshold": ami_max,
            "clusters": clusters,
            "conditions": conditions,
        }

        reasoning_lines: List[str] = []
        reasoning_lines.append(f"Adjusted Mutual Information: {ami_value:.4f}")
        threshold_desc = []
        if ami_min is not None:
            threshold_desc.append(f"≥{ami_min:.4f}")
        if ami_max is not None:
            threshold_desc.append(f"≤{ami_max:.4f}")
        if threshold_desc:
            reasoning_lines.append(f"Threshold: {' and '.join(threshold_desc)}")
        reasoning_lines.append(f"Result: {'PASS' if passed else 'FAIL'}")

        return GraderResult(
            passed=passed,
            metrics=metrics,
            reasoning="\n".join(reasoning_lines),
            agent_answer=agent_answer,
        )
