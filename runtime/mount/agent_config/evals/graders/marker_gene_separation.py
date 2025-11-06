from .base import BinaryGrader, GraderResult
from eval_types import TestResult

class MarkerGeneSeparationGrader(BinaryGrader):
    def evaluate(self, test_result: TestResult, config: dict) -> GraderResult:
        scoring = config.get("scoring", {})
        thresholds = scoring.get("pass_thresholds", {})
        mean_auroc_threshold = thresholds.get("mean_auroc", 0.85)
        fraction_high_threshold = thresholds.get("fraction_high", 0.70)
        per_gene_cutoff = thresholds.get("per_gene_cutoff", 0.80)

        agent_answer = self.extract_answer_from_tags(test_result.conversation_history)

        if agent_answer is None:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Agent did not provide answer in <EVAL_ANSWER> tags",
                agent_answer=None
            )

        if "per_gene_stats" not in agent_answer:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Agent answer missing required field: per_gene_stats",
                agent_answer=agent_answer
            )

        if "mean_auroc" not in agent_answer:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="Agent answer missing required field: mean_auroc",
                agent_answer=agent_answer
            )

        per_gene_stats = agent_answer["per_gene_stats"]
        agent_mean_auroc = agent_answer["mean_auroc"]

        if not isinstance(per_gene_stats, list):
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="per_gene_stats must be a list",
                agent_answer=agent_answer
            )

        num_genes = len(per_gene_stats)
        if num_genes == 0:
            return GraderResult(
                passed=False,
                metrics={},
                reasoning="per_gene_stats is empty",
                agent_answer=agent_answer
            )

        gene_aurocs = {}
        for stat in per_gene_stats:
            if not isinstance(stat, dict):
                return GraderResult(
                    passed=False,
                    metrics={},
                    reasoning="Each element in per_gene_stats must be a dict with 'gene' and 'auroc'",
                    agent_answer=agent_answer
                )
            if "gene" not in stat or "auroc" not in stat:
                return GraderResult(
                    passed=False,
                    metrics={},
                    reasoning="Each element in per_gene_stats must have 'gene' and 'auroc' fields",
                    agent_answer=agent_answer
                )
            gene_aurocs[stat["gene"]] = stat["auroc"]

        computed_mean_auroc = sum(gene_aurocs.values()) / len(gene_aurocs)

        high_auroc_genes = [gene for gene, auroc in gene_aurocs.items() if auroc >= per_gene_cutoff]
        low_auroc_genes = [gene for gene, auroc in gene_aurocs.items() if auroc < per_gene_cutoff]
        fraction_high = len(high_auroc_genes) / num_genes

        mean_auroc_pass = agent_mean_auroc >= mean_auroc_threshold
        fraction_high_pass = fraction_high >= fraction_high_threshold

        passed = mean_auroc_pass and fraction_high_pass

        metrics = {
            "num_genes": num_genes,
            "mean_auroc_agent": agent_mean_auroc,
            "mean_auroc_computed": computed_mean_auroc,
            "mean_auroc_threshold": mean_auroc_threshold,
            "fraction_high": fraction_high,
            "fraction_high_threshold": fraction_high_threshold,
            "per_gene_cutoff": per_gene_cutoff,
            "num_high_auroc_genes": len(high_auroc_genes),
            "num_low_auroc_genes": len(low_auroc_genes),
            "high_auroc_genes": sorted(high_auroc_genes),
            "low_auroc_genes": sorted(low_auroc_genes),
            "mean_auroc_pass": mean_auroc_pass,
            "fraction_high_pass": fraction_high_pass,
            "per_gene_aurocs": gene_aurocs,
        }

        reasoning = self._format_separation_reasoning(
            num_genes,
            agent_mean_auroc,
            computed_mean_auroc,
            mean_auroc_threshold,
            fraction_high,
            fraction_high_threshold,
            per_gene_cutoff,
            gene_aurocs,
            high_auroc_genes,
            low_auroc_genes,
            mean_auroc_pass,
            fraction_high_pass,
            passed
        )

        return GraderResult(
            passed=passed,
            metrics=metrics,
            reasoning=reasoning,
            agent_answer=agent_answer
        )

    def _format_separation_reasoning(self, num_genes, agent_mean, computed_mean,
                                     mean_threshold, fraction_high, fraction_threshold,
                                     per_gene_cutoff, gene_aurocs, high_genes, low_genes,
                                     mean_pass, fraction_pass, passed):
        lines = []
        lines.append(f"Marker Gene Separation: {'PASS' if passed else 'FAIL'}")
        lines.append("")

        mean_check = "✓" if mean_pass else "✗"
        lines.append(f"Mean AUROC: {agent_mean:.3f} {mean_check} (threshold: {mean_threshold:.3f})")

        fraction_check = "✓" if fraction_pass else "✗"
        lines.append(f"Fraction High AUROC (≥{per_gene_cutoff:.2f}): {fraction_high:.3f} ({len(high_genes)}/{num_genes}) {fraction_check} (threshold: {fraction_threshold:.3f})")

        lines.append("")
        lines.append("Per-gene AUROC scores:")
        for gene in sorted(gene_aurocs.keys(), key=lambda g: gene_aurocs[g], reverse=True):
            auroc = gene_aurocs[gene]
            check = "✓" if auroc >= per_gene_cutoff else "✗"
            lines.append(f"  {check} {gene}: {auroc:.3f}")

        if abs(agent_mean - computed_mean) > 0.001:
            lines.append("")
            lines.append(f"Note: Agent reported mean {agent_mean:.3f}, computed mean is {computed_mean:.3f}")

        lines.append("")
        lines.append(f"Result: {'PASS' if passed else 'FAIL'}")
        if not passed:
            failures = []
            if not mean_pass:
                failures.append(f"Mean AUROC {agent_mean:.3f} < {mean_threshold:.3f}")
            if not fraction_pass:
                failures.append(f"Fraction high {fraction_high:.3f} < {fraction_threshold:.3f}")
            lines.append(f"Reasons: {'; '.join(failures)}")

        return "\n".join(lines)
