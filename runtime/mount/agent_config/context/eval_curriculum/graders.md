# Grader Configs

## numeric_tolerance

QC metrics, counts, normalized values.

```json
{
  "grader": {
    "type": "numeric_tolerance",
    "config": {
      "ground_truth": {"cells_after_filtering": 1374915},
      "tolerances": {
        "cells_after_filtering": {"type": "absolute", "value": 50}
      }
    }
  }
}
```

Tolerance types: `absolute` (±50) or `relative` (±0.05 = 5%).

**Calibration:**
- Tight (±0.004%) when answer is stable (clean data, low-quality tail removal)
- Wide (±20%) when multiple valid methods exist

## multiple_choice

Interpretation, pattern identification.

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

**Rules:**
- All distractors must be biologically plausible
- Remove descriptive labels ("Cluster 3" not "Cluster 3 (bone formation)")
- Question must require data inspection

## distribution_comparison

Cell type proportions, population compositions.

```json
{
  "grader": {
    "type": "distribution_comparison",
    "config": {
      "ground_truth": {
        "cell_type_distribution": {"TAL": 14.81, "Fib": 11.66, "PTS1": 8.89}
      },
      "tolerances": {
        "cell_type_percentages": {"type": "absolute", "value": 5.0}
      }
    }
  }
}
```

Typically ±3-5 percentage points per cell type.

## marker_gene_precision_recall

Marker discovery, DE results, top-K gene lists.

```json
{
  "grader": {
    "type": "marker_gene_precision_recall",
    "config": {
      "canonical_markers": ["COL1A1", "COL1A2", "SPP1", "SPARC"],
      "scoring": {
        "pass_thresholds": {
          "recall_at_k": 0.5,
          "precision_at_k": 0
        }
      }
    }
  }
}
```

Recall threshold = fraction of canonical markers that must appear in agent's list.

## label_set_jaccard

Multi-select, set matching, pathway selection.

```json
{
  "grader": {
    "type": "label_set_jaccard",
    "config": {
      "ground_truth": ["A", "C", "E"],
      "threshold": 0.67
    }
  }
}
```

Jaccard = intersection / union. Penalizes both missing and extra items.

## Decision Tree

```
Single number? → numeric_tolerance
Single choice among options? → multiple_choice
Multiple selections? → label_set_jaccard
Gene list (ranked)? → marker_gene_precision_recall
Gene list (unordered)? → label_set_jaccard
Cell type proportions? → distribution_comparison
```
