import numpy as np
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests

from .base import BinaryGrader, GraderResult
from eval_types import TestResult

class ProportionConsistencyGrader(BinaryGrader):
    def evaluate(self, test_result: TestResult, config: dict) -> GraderResult:
        test_config = config.get("test", {})
        scoring = config.get("scoring", {})
        thresholds = scoring.get("pass_thresholds", {})
        chi2_q_min = thresholds.get("chi2_q_min", 0.05)
        cramers_v_max = thresholds.get("cramers_v_max", 0.15)

        agent_answer = self.extract_answer_from_tags(test_result.conversation_history)

        if agent_answer is None:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Agent did not provide answer in <EVAL_ANSWER> tags",
                agent_answer=None
            )

        if "cell_types" not in agent_answer:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Agent answer missing required field: cell_types",
                agent_answer=agent_answer
            )

        if "samples" not in agent_answer:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Agent answer missing required field: samples",
                agent_answer=agent_answer
            )

        cell_types = agent_answer["cell_types"]
        samples_data = agent_answer["samples"]

        if not isinstance(cell_types, list) or len(cell_types) == 0:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="cell_types must be a non-empty list",
                agent_answer=agent_answer
            )

        if not isinstance(samples_data, list) or len(samples_data) == 0:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="samples must be a non-empty list",
                agent_answer=agent_answer
            )

        contingency_table = []
        sample_names = []

        for sample_data in samples_data:
            if not isinstance(sample_data, dict):
                return GraderResult(
                    passed=False,
                    metrics={},
                    reasoning="Each sample entry must be a dict with 'sample' and 'proportions'",
                    agent_answer=agent_answer
                )

            if "sample" not in sample_data or "proportions" not in sample_data:
                return GraderResult(
                    passed=False,
                    metrics={},
                    reasoning="Each sample entry must have 'sample' and 'proportions' fields",
                    agent_answer=agent_answer
                )

            sample_name = sample_data["sample"]
            proportions = sample_data["proportions"]
            sample_names.append(sample_name)

            row = [proportions.get(ct, 0.0) for ct in cell_types]
            contingency_table.append(row)

        contingency_array = np.array(contingency_table)

        if contingency_array.shape[0] < 2:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Need at least 2 samples to perform chi-square test",
                agent_answer=agent_answer
            )

        try:
            chi2_stat, p_value, dof, expected = chi2_contingency(contingency_array)
        except Exception as e:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning=f"Failed to compute chi-square test: {e}",
                agent_answer=agent_answer
            )

        rejected, q_values, _, _ = multipletests([p_value], method='fdr_bh')
        q_value = q_values[0]

        n = contingency_array.sum()
        min_dim = min(contingency_array.shape[0] - 1, contingency_array.shape[1] - 1)
        cramers_v = np.sqrt(chi2_stat / (n * min_dim)) if min_dim > 0 and n > 0 else 0.0

        q_pass = q_value >= chi2_q_min
        cramers_pass = cramers_v <= cramers_v_max

        passed = q_pass and cramers_pass

        metrics = {
            "num_samples": len(sample_names),
            "num_cell_types": len(cell_types),
            "chi2_statistic": chi2_stat,
            "p_value": p_value,
            "q_value": q_value,
            "degrees_of_freedom": dof,
            "cramers_v": cramers_v,
            "chi2_q_min_threshold": chi2_q_min,
            "cramers_v_max_threshold": cramers_v_max,
            "q_pass": q_pass,
            "cramers_pass": cramers_pass,
            "cell_types": cell_types,
            "sample_names": sample_names,
        }

        reasoning = self._format_proportion_consistency_reasoning(
            len(sample_names),
            len(cell_types),
            chi2_stat,
            p_value,
            q_value,
            dof,
            cramers_v,
            chi2_q_min,
            cramers_v_max,
            q_pass,
            cramers_pass,
            passed,
            sample_names,
            cell_types
        )

        return GraderResult(
            passed=passed,
            metrics=metrics,
            reasoning=reasoning,
            agent_answer=agent_answer
        )

    def _format_proportion_consistency_reasoning(self, num_samples, num_cell_types,
                                                 chi2_stat, p_value, q_value, dof,
                                                 cramers_v, q_min, v_max,
                                                 q_pass, cramers_pass, passed,
                                                 sample_names, cell_types):
        lines = []
        lines.append(f"Proportion Consistency (Chi-Square Test): {'PASS' if passed else 'FAIL'}")
        lines.append("")

        lines.append(f"Samples analyzed: {num_samples} ({', '.join(sample_names)})")
        lines.append(f"Cell types: {num_cell_types} ({', '.join(cell_types)})")
        lines.append("")

        lines.append("Chi-Square Test of Independence:")
        lines.append(f"  χ² statistic: {chi2_stat:.4f}")
        lines.append(f"  Degrees of freedom: {dof}")
        lines.append(f"  p-value: {p_value:.6f}")

        q_check = "✓" if q_pass else "✗"
        lines.append(f"  q-value (Benjamini-Hochberg): {q_value:.6f} {q_check} (threshold: ≥{q_min})")

        lines.append("")
        cramers_check = "✓" if cramers_pass else "✗"
        lines.append(f"Effect Size (Cramér's V): {cramers_v:.4f} {cramers_check} (threshold: ≤{v_max})")

        lines.append("")
        lines.append(f"Result: {'PASS' if passed else 'FAIL'}")
        if not passed:
            failures = []
            if not q_pass:
                failures.append(f"q-value {q_value:.6f} < {q_min} (significant heterogeneity detected)")
            if not cramers_pass:
                failures.append(f"Cramér's V {cramers_v:.4f} > {v_max} (effect size too large)")
            lines.append(f"Reasons: {'; '.join(failures)}")

        return "\n".join(lines)
