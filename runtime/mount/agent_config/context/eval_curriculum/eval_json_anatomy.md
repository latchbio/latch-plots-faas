# Eval JSON Anatomy

This document describes the structure of an eval JSON file, with emphasis on required and recommended metadata fields.

## Complete Structure

```json
{
  "id": "string (required)",
  "task": "string (required)",
  "data_node": "string | string[] (required)",
  "grader": {
    "type": "string (required)",
    "config": { ... }
  },
  "notes": "string (optional)",
  "metadata": {
    "task": "string (required)",
    "time_horizon": "string (required)",
    "kit": "string (required)",
    "eval_type": "string (required)"
  }
}
```

## Field Reference

### Core Fields

#### `id` (required)
Unique identifier for the eval. Should be descriptive and use snake_case.

```json
"id": "xenium_qc_filter_min_umi_counts"
```

**Conventions:**
- Start with platform/kit name when relevant
- Include task category
- Be specific about what's being tested

#### `task` (required)
Natural language description of what the agent should do. Must include:
- Clear objective
- Expected output format with field names
- Any constraints or requirements

```json
"task": "Filter cells with very low total UMI counts. Return a JSON object with field: cells_after_filtering (number remaining after filter)."
```

**Best practices:**
- Explicitly specify output field names
- Use `Return EXACTLY:` before the `<EVAL_ANSWER>` block (required by linter)
- For numeric fields, use unquoted placeholders: `{"count": <int>}` not `{"count": "<value>"}`
- Avoid ambiguity in what constitutes success

**Answer format convention:**
```
Return EXACTLY:

<EVAL_ANSWER>
{"field_name": <int>}
</EVAL_ANSWER>
```

**Multiple choice format:**
```
Return EXACTLY:

<EVAL_ANSWER>
{"answer": "<letter>"}
</EVAL_ANSWER>

where <letter> is A, B, C, or D.
```
Use `<letter>` as the placeholder for multiple choice answers (not `<value>` or `X`).

The agent must include the `<EVAL_ANSWER>` tags in their response. The harness extracts JSON from within these tags.

#### `data_node` (required)
Latch data node URL(s) pointing to the input data.

```json
"data_node": "latch://12345.node"
```

Or for multiple files:
```json
"data_node": [
  "latch://12345.node",
  "latch://67890.node"
]
```

Use stable data ID format (`latch://12345.node`), not path-based URLs.

#### `grader` (required)
Specifies the grader type and configuration.

```json
"grader": {
  "type": "numeric_tolerance",
  "config": {
    "ground_truth": { "cells_after_filtering": 1374915 },
    "tolerances": {
      "cells_after_filtering": { "type": "absolute", "value": 50 }
    }
  }
}
```

See [Graders](#graders) section in README for all grader types.

#### `notes` (optional but recommended)
Internal notes not shown to the agent. Use this field to capture:
- How to pass the test (the expected solution approach)
- Tolerance calibration rationale
- Design decisions and known edge cases

Since evals should be tests without instructions, document the solution path here. Additional markdown files alongside the eval JSON are encouraged for detailed observations or methodology.

```json
"notes": "Tolerance is tight because this dataset looks clean and the paper uses very low QC cutoffs (<10 UMIs). In a local test, the agent used an overly high cutoff (100 UMIs) without checking the UMI distribution."
```

### Metadata Fields

The `metadata` object captures information for filtering, analysis, and benchmark organization.

#### `metadata.task` (required)
Task category. Standard values:
- `qc` - Quality control filtering
- `normalization` - Data normalization/transformation
- `clustering` - Clustering analysis
- `cell_typing` - Cell type annotation/classification
- `differential_expression` - DE analysis
- `dimension_reduction` - PCA, UMAP, etc.
- `spatial_analysis` - Spatial statistics, neighborhood analysis

```json
"metadata": {
  "task": "qc"
}
```

#### `metadata.time_horizon` (required)
Semantic grouping for expected computational complexity. Used for filtering and analysis, **not** for setting actual timeouts.

- `small` - Simple lookups, basic statistics, small data operations
- `medium` - Standard analysis workflows (PCA, clustering, DE)
- `large` - Complex pipelines, large datasets, multiple analysis steps

```json
"metadata": {
  "time_horizon": "small"
}
```

**Note:** This is a semantic label for grouping evals by complexity class. To set actual execution timeouts, use `timeout_s`.

#### `metadata.kit` (required)
Spatial transcriptomics platform:
- `xenium` - 10x Genomics Xenium
- `visium` - 10x Genomics Visium
- `merfish` / `vizgen` - Vizgen MERFISH
- `cosmx` - NanoString CosMx
- `seeker` / `takara` - Takara PLL-Seq / Seeker
- `atlasxomics` - AtlasXomics DBiT-seq
- `curio` - Curio Seeker

```json
"metadata": {
  "kit": "xenium"
}
```

#### `metadata.eval_type` (required)
Categorizes the eval by what it tests and how durability criteria are applied:
- `scientific` - Agent selects appropriate method for a scientific goal (wide durability tolerance)
- `procedural` - Agent applies a named method correctly (moderate durability tolerance)
- `observational` - Concept exposure or training data generation (relaxed durability)

```json
"metadata": {
  "eval_type": "scientific"
}
```

#### `metadata.timeout_s` (optional)
Override the default eval timeout (600 seconds). Use for evals that legitimately require more time due to large datasets or complex computations.

```json
"metadata": {
  "timeout_s": 900
}
```

**When to use:**
- Eval consistently times out despite productive agent work
- Large dataset requires extended download or processing time
- Multi-step analysis that legitimately takes longer

**Finding the right timeout:**
1. Run the eval with increasing `default_timeout_s` values in the workflow (e.g., 600, 900, 1200)
2. Check the trajectory to see how long successful runs take
3. Set `timeout_s` to ~1.5x the typical successful duration
4. If an eval needs >15 minutes, consider whether the task scope is too broad

**Important:** Always verify the task is well-scoped before increasing timeout. A timeout often indicates the agent is stuck or the task is ambiguous, not that more time is needed.

**Choosing the right type:**

| If the task... | Use |
|----------------|-----|
| Describes a scientific goal without specifying method | `scientific` |
| Names a specific method but not parameters | `procedural` |
| Is for concept exposure or training data | `observational` |

**Example distinction:**
- `scientific`: "Filter low-quality cells" (agent picks method + params)
- `procedural`: "Apply negative probe filtering" (method named, params not specified)
- `observational`: "Explore the negative probe background signal" (loose testing)

## Complete Example

```json
{
  "id": "xenium_kidney_frpt_vs_healthy_pt_deg",
  "task": "Identify the top 10 differentially expressed genes between FR-PT (failed-repair proximal tubule) and healthy PT cells. Use rank_genes_groups with Wilcoxon test. Return EXACTLY:\n<EVAL_ANSWER>\n{\"top_marker_genes\": [\"GENE1\", \"GENE2\", ...]}\n</EVAL_ANSWER>",
  "data_node": "latch://38438.account/Scratch/xenium/evals/xenium_mouse_kidney_annotated.h5ad",
  "grader": {
    "type": "marker_gene_precision_recall",
    "config": {
      "canonical_markers": ["Havcr1", "Vcam1", "Krt20", "Dcdc2a", "Ccl2"],
      "scoring": {
        "pass_thresholds": {
          "precision_at_k": 0.0,
          "recall_at_k": 0.6
        }
      }
    }
  },
  "notes": "Canonical markers from Lake et al. 2023 kidney injury atlas. Recall threshold set at 60% (3/5 markers) to allow for dataset-specific variation.",
  "metadata": {
    "task": "differential_expression",
    "time_horizon": "small",
    "kit": "xenium",
    "eval_type": "procedural"
  }
}
```

## Validation Checklist

Before submitting an eval, verify:

- [ ] `id` is unique and descriptive
- [ ] `task` explicitly specifies output format and field names
- [ ] `data_node` uses stable ID format (not path-based)
- [ ] `grader.type` matches expected answer structure
- [ ] `grader.config` has calibrated tolerances with rationale
- [ ] `metadata.task` uses standard category
- [ ] `metadata.time_horizon` is set
- [ ] `metadata.kit` matches the platform
- [ ] `metadata.eval_type` is set (scientific, procedural, or observational)
- [ ] `notes` explains any non-obvious design decisions
