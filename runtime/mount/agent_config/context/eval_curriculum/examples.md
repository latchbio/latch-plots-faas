# Eval Examples

## Good: Xenium QC - UMI Filtering

```json
{
  "id": "xenium_qc_filter_min_umi_counts",
  "task": "Filter cells with very low total UMI counts. Return JSON: {cells_after_filtering: <number>}",
  "data_node": "latch://38438.account/Scratch/xenium/Output/xenium_combined_12samples.h5ad",
  "grader": {
    "type": "numeric_tolerance",
    "config": {
      "ground_truth": {"cells_after_filtering": 1374915},
      "tolerances": {"cells_after_filtering": {"type": "absolute", "value": 50}}
    }
  },
  "notes": "Tolerance tight because dataset is clean. Agent must inspect UMI distribution - naive cutoffs (100 UMIs) filtered 40% of cells.",
  "metadata": {"task": "qc", "kit": "xenium"}
}
```

**Why good:**
- Verifiable: clear numeric output, deterministic grader
- Durable: clean data means any reasonable threshold (5-15 UMIs) gives same result
- Anti-shortcut: must inspect distribution; naive cutoffs fail

## Good: Kidney Cell Type Distribution

```json
{
  "id": "kidney_major_20populations_distribution",
  "task": "Assign cells to 20 populations, return JSON: {total_cells: <n>, cell_type_distribution: {TypeA: <pct>, ...}}",
  "grader": {
    "type": "distribution_comparison",
    "config": {
      "ground_truth": {"cell_type_distribution": {"TAL": 14.81, "Fib": 11.66, ...}},
      "tolerances": {"cell_type_percentages": {"type": "absolute", "value": 3.0}}
    }
  }
}
```

**Why good:** Tests biology (cell composition), not implementation details.

## Bad: Subjective Task

```json
{
  "task": "Find the most interesting differentially expressed genes"
}
```

**Problem:** "Interesting" has no quantitative definition. Not verifiable.

**Fix:** "Return top 10 DE genes by log fold change as JSON: {top_genes: [...]}"

## Bad: Hyperparameter Sensitive

```json
{
  "task": "How many clusters does Leiden produce?",
  "ground_truth": {"n_clusters": 12}
}
```

**Problem:** Cluster count depends on resolution and random seed. Not durable.

**Fix:** Test interpretation ("Which cluster shows bone formation markers?") instead of count.

## Bad: Prior Knowledge Solvable

```json
{
  "task": "What marker gene distinguishes podocytes from other kidney cells?",
  "options": ["NPHS2", "COL1A1", "AQP2", "SLC12A1"]
}
```

**Problem:** Answerable from textbook knowledge without data. Not anti-shortcut.

**Fix:** "Which gene shows highest AUROC for podocyte vs rest in this dataset?"

## Bad: Answer in Question

```json
{
  "task": "Which cluster contains bone formation cells?",
  "options": ["Cluster 1 (immune)", "Cluster 2 (bone formation)", "Cluster 3 (epithelial)"]
}
```

**Problem:** Option label leaks answer. Not anti-shortcut.

**Fix:** Remove descriptive labels: "Cluster 1", "Cluster 2", "Cluster 3"
