# Eval Design Curriculum

<detection>
Use when user wants to:
- Create evaluation problems for SpatialBench
- Design graders or test cases
- Review/improve existing evals
</detection>

## Three Design Principles

Every eval must satisfy:

### 1. Verifiability
- Deterministic grader with unambiguous pass/fail
- Output format explicitly specified with JSON field names
- No subjective interpretation ("interesting genes" = bad)

### 2. Scientific Durability
- Answer stable across random seeds, preprocessing choices, library versions
- Tests biology, not implementation artifacts
- Multiple valid methods converge on same answer
- Tolerance calibrated to accept valid variation

### 3. Anti-Shortcut
- Cannot be solved by prior knowledge alone
- Requires actual data inspection/computation
- Answer not leaked in question wording
- Distractors (if multiple choice) require data analysis to eliminate

## Eval Types

| Type | Task Style | Durability |
|------|-----------|------------|
| **Scientific** | "Apply comprehensive QC" - agent picks methods AND parameters | Wide tolerance |
| **Procedural** | "Apply negative probe filtering" - method specified, agent picks parameters | Moderate tolerance |
| **Observational** | Exposure/training data generation | Relaxed |

## JSON Structure

```json
{
  "id": "xenium_qc_filter_min_umi_counts",
  "task": "Filter cells with low UMI counts. Return JSON: {cells_after_filtering: <number>}",
  "data_node": "latch://...",
  "grader": {
    "type": "numeric_tolerance",
    "config": { ... }
  },
  "notes": "Document curator reasoning, observed failure modes",
  "metadata": {
    "task": "qc",
    "time_horizon": "small",
    "kit": "xenium"
  }
}
```

## Grader Types

See `graders.md` for detailed configs.

| Grader | Use Case | Answer Format |
|--------|----------|---------------|
| `numeric_tolerance` | Counts, metrics, ratios | `{"field": <number>}` |
| `multiple_choice` | Interpretation, patterns | `{"answer": "A"}` |
| `distribution_comparison` | Cell type proportions | `{"distribution": {"TypeA": <pct>, ...}}` |
| `marker_gene_precision_recall` | Gene lists, DE results | `{"top_marker_genes": [...]}` |
| `label_set_jaccard` | Multi-select, set matching | `{"labels": ["A", "C"]}` |

## Common Mistakes

| Mistake | Principle Violated | Fix |
|---------|-------------------|-----|
| "Find interesting..." | Verifiability | Define quantitative criteria |
| No output format | Verifiability | Specify JSON fields in task |
| Testing cluster count | Durability | Test interpretation instead |
| Testing coordinates | Durability | Test relative positions |
| Textbook question | Anti-Shortcut | Require dataset-specific computation |
