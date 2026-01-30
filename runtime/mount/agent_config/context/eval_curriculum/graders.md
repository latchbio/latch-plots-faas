# Grader Types

## `numeric_tolerance`

**Use case:** QC metrics, counts, normalized values, ratios

**Answer format:** `{"field_name": <number>}`

**How it works:** Checks if the agent's numeric answer is within a specified tolerance of ground truth.

**Tolerance types:**
- `absolute` - value must be within ±N of expected (default)
- `relative` - value must be within ±N% of expected
- `min` - value must be >= threshold
- `max` - value must be <= threshold
- Asymmetric absolute tolerances using `lower` and `upper` bounds

**Example config:**
```json
{
  "grader": {
    "type": "numeric_tolerance",
    "config": {
      "ground_truth": { "cells_after_filtering": 1374915 },
      "tolerances": {
        "cells_after_filtering": { "type": "absolute", "value": 50 }
      }
    }
  }
}
```

**Example with min/max thresholds:**
```json
{
  "grader": {
    "type": "numeric_tolerance",
    "config": {
      "ground_truth": { "score": 0.85 },
      "tolerances": {
        "score": { "type": "min", "value": 0.80 }
      }
    }
  }
}
```

**Example with asymmetric tolerance:**
```json
{
  "grader": {
    "type": "numeric_tolerance",
    "config": {
      "ground_truth": { "count": 100 },
      "tolerances": {
        "count": { "type": "absolute", "lower": 10, "upper": 20 }
      }
    }
  }
}
```

**When to use:** When the answer is a single number and you can define acceptable variation.

---

## `multiple_choice`

**Use case:** Interpretation questions, pattern identification, qualitative comparisons

**Answer format:** `{"answer": "A"}` (case-insensitive)

**How it works:** Exact match against the correct answer letter. Supports multiple correct answers.

**Example config:**
```json
{
  "grader": {
    "type": "multiple_choice",
    "config": {
      "correct_answer": "B"
    }
  }
}
```

**Example with multiple correct answers:**
```json
{
  "grader": {
    "type": "multiple_choice",
    "config": {
      "correct_answers": ["B", "C"]
    }
  }
}
```

**When to use:** When testing interpretation or understanding rather than computation. All distractors must be biologically plausible.

---

## `distribution_comparison`

**Use case:** Cell type proportions, population compositions

**Answer format:** `{"cell_type_distribution": {"TypeA": <pct>, "TypeB": <pct>, ...}}`

**How it works:** Checks each value in the distribution against per-field tolerances. All fields must be within tolerance to pass.

**Example config:**
```json
{
  "grader": {
    "type": "distribution_comparison",
    "config": {
      "ground_truth": {
        "cell_type_distribution": {
          "TAL": 14.81,
          "Fib": 11.66,
          "PTS1": 8.89
        }
      },
      "tolerances": {
        "cell_type_percentages": { "type": "absolute", "value": 5.0 }
      }
    }
  }
}
```

**When to use:** When testing multi-class annotation where each category has an expected proportion.

---

## `marker_gene_precision_recall`

**Use case:** Marker gene discovery, differential expression, top-K gene lists

**Answer format:** `{"top_marker_genes": ["Gene1", "Gene2", ...]}` (configurable via `answer_field`)

**How it works:** Compares agent's gene list against a canonical marker set. Computes recall (fraction of canonical markers recovered) and precision. Supports both flat lists and per-celltype dictionaries.

**Defaults:** `precision_at_k: 0.60`, `recall_at_k: 0.50`

**Example config:**
```json
{
  "grader": {
    "type": "marker_gene_precision_recall",
    "config": {
      "canonical_markers": ["COL1A1", "COL1A2", "SPP1", "SPARC", "BGLAP", "IBSP"],
      "answer_field": "top_marker_genes",
      "scoring": {
        "pass_thresholds": {
          "recall_at_k": 0.50,
          "precision_at_k": 0.60
        }
      }
    }
  }
}
```

**Example config for per-celltype evaluation:**
```json
{
  "grader": {
    "type": "marker_gene_precision_recall",
    "config": {
      "canonical_markers": {
        "T_cells": ["CD3D", "CD3E", "CD4"],
        "B_cells": ["CD19", "MS4A1", "CD79A"]
      },
      "scoring": {
        "pass_thresholds": {
          "min_recall_per_celltype": 0.50,
          "min_celltypes_passing": 2
        }
      }
    }
  }
}
```

**When to use:** For discovery tasks where you have a validated set of expected markers. Recall threshold determines how many canonical markers must appear in the agent's list.

---

## `label_set_jaccard`

**Aliases:** `jaccard_label_set`

**Use case:** Multi-select questions, set matching, pathway selection

**Answer format:** `{"cell_types_predicted": ["A", "C", "D"]}` (configurable via `answer_field`)

**How it works:** Computes Jaccard index (intersection / union) between agent's set and ground truth set. Threshold determines minimum similarity to pass.

**Default:** `pass_threshold: 0.90`, `answer_field: "cell_types_predicted"`

**Example config:**
```json
{
  "grader": {
    "type": "label_set_jaccard",
    "config": {
      "ground_truth_labels": ["A", "C", "E"],
      "answer_field": "cell_types_predicted",
      "scoring": {
        "pass_threshold": 0.67
      }
    }
  }
}
```

**When to use:** When multiple answers are correct and order doesn't matter. Penalizes both missing items and extra incorrect items.

---

## `marker_gene_separation`

**Use case:** Evaluating how well marker genes discriminate their target cell type from others

**Answer format:**
```json
{
  "mean_auroc": 0.92,
  "per_gene_stats": [
    {"gene": "CD3D", "auroc": 0.95},
    {"gene": "CD3E", "auroc": 0.89}
  ]
}
```

**How it works:** Evaluates marker gene quality by checking mean AUROC and the fraction of genes exceeding a per-gene AUROC cutoff.

**Defaults:** `mean_auroc: 0.85`, `fraction_high: 0.70`, `per_gene_cutoff: 0.80`

**Example config:**
```json
{
  "grader": {
    "type": "marker_gene_separation",
    "config": {
      "scoring": {
        "pass_thresholds": {
          "mean_auroc": 0.85,
          "fraction_high": 0.70,
          "per_gene_cutoff": 0.80
        }
      }
    }
  }
}
```

**When to use:** When you need to verify that identified markers actually discriminate the target population well.

---

## `spatial_adjacency`

**Use case:** Spatial transcriptomics, cell-cell proximity analysis

**Answer format:**
```json
{
  "median_ic_to_pc_um": 12.5,
  "p90_ic_to_pc_um": 45.0,
  "pct_ic_within_15um": 72.0,
  "pct_ic_mixed_within_55um": 85.0,
  "adjacency_pass": true
}
```

**How it works:** Evaluates spatial relationships between cell types (e.g., immune cells to parenchymal cells) using distance metrics and proximity thresholds.

**Defaults:** `max_median_ic_to_pc_um: 25.0`, `max_p90_ic_to_pc_um: 80.0`, `min_pct_ic_within_15um: 60.0`, `min_pct_ic_mixed_within_55um: 60.0`

**Example config:**
```json
{
  "grader": {
    "type": "spatial_adjacency",
    "config": {
      "scoring": {
        "pass_thresholds": {
          "max_median_ic_to_pc_um": 25.0,
          "max_p90_ic_to_pc_um": 80.0,
          "min_pct_ic_within_15um": 60.0,
          "min_pct_ic_mixed_within_55um": 60.0
        }
      }
    }
  }
}
```

**When to use:** For spatial transcriptomics tasks where cell-cell proximity or neighborhood composition matters.
