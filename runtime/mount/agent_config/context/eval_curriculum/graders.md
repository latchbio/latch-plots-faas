## Grader Types

### 2.1 `numeric_tolerance`

**Use case:** QC metrics, counts, normalized values, ratios

**Answer format:** `{"field_name": <number>}`

**How it works:** Checks if the agent's numeric answer is within a specified tolerance of ground truth. Tolerance can be absolute (±50) or relative (±5%).

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

**When to use:** When the answer is a single number and you can define acceptable variation. See Examples 1, 2, 3, 10.

---

### 2.2 `multiple_choice`

**Use case:** Interpretation questions, pattern identification, qualitative comparisons

**Answer format:** `{"answer": "A"}` (case-insensitive)

**How it works:** Exact match against the correct answer letter.

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

**When to use:** When testing interpretation or understanding rather than computation. All distractors must be biologically plausible. See Examples 4, 9.

---

### 2.3 `distribution_comparison`

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

**When to use:** When testing multi-class annotation where each category has an expected proportion. See Example 6.

---

### 2.4 `marker_gene_precision_recall`

**Use case:** Marker gene discovery, differential expression, top-K gene lists

**Answer format:** `{"top_marker_genes": ["Gene1", "Gene2", ...]}`

**How it works:** Compares agent's gene list against a canonical marker set. Computes recall (fraction of canonical markers recovered) and optionally precision.

**Example config:**
```json
{
  "grader": {
    "type": "marker_gene_precision_recall",
    "config": {
      "canonical_markers": ["COL1A1", "COL1A2", "SPP1", "SPARC", "BGLAP", "IBSP"],
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

**When to use:** For discovery tasks where you have a validated set of expected markers. Recall threshold determines how many canonical markers must appear in the agent's list. See Examples 5, 7, 8.

---

### 2.5 `label_set_jaccard`

**Use case:** Multi-select questions, set matching, pathway selection

**Answer format:** `{"labels": ["A", "C", "D"]}` or `{"pathways": ["apoptosis", "inflammation"]}`

**How it works:** Computes Jaccard index (intersection / union) between agent's set and ground truth set. Threshold determines minimum similarity to pass.

**Example config:**
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

**When to use:** When multiple answers are correct and order doesn't matter. Penalizes both missing items and extra incorrect items.

