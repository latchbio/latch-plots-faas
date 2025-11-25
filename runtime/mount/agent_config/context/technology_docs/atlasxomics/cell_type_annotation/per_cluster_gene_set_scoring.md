## Cluster-Based Gene Set Scoring - Concise Guide

### Method Overview

Identify cell types using fold change enrichment scoring: `FC = mean_in_cluster / mean_not_in_cluster`

**Key principle:** FC > 1 = enriched (assign), FC ≤ 1 = not enriched (mark Unknown)

---

## Workflow

1. **Define cell types** — Collect 10+ marker genes per type from CellGuide or literature
   - **CRITICAL:** Never omit expected cell types. If CellGuide lacks markers, supplement from literature

2. **Check markers** — Verify ≥5 markers per cell type exist in `adata.var_names`

3. **Compute fold change** — For each cluster × cell type:
   ```python
   FC = mean_expression_in_cluster / mean_expression_in_non_cluster_cells
   ```

4. **Assign cell types:**
   - If max FC > threshold (default 1.01): assign that cell type
   - If max FC ≤ threshold: mark as "Unknown"
   - Store in `adata.obs['cell_type']`

5. **⚠️ CRITICAL: Evaluate your work** — Read `technology_docs/atlasxomics/cell_type_annotation/evals.md`
   - Compute ALL required metrics: proportions, purity, spatial coherence, marker enrichment, confidence, sample consistency, condition effects

6. **Revise your work and add corrective steps** if results are weak

---

## Critical Parameters

**Fold Change Threshold:**
- **1.01** (>1% enrichment) — **recommended default**
- 1.0 — permissive (any enrichment)
- 1.05 — moderate (>5% enrichment)
- 1.1 — stringent (>10% enrichment)

**Confidence:**
- High: Top FC - second FC > 0.1
- Moderate: Separation 0.05-0.1
- Low/Unknown: Top FC ≤ threshold or separation < 0.05

---

## Implementation

```python
import numpy as np
import pandas as pd

clusters = adata.obs['cluster'].unique()
fold_change_scores = pd.DataFrame(index=clusters, columns=marker_dict.keys())

for cluster in clusters:
    cluster_mask = adata.obs['cluster'] == cluster
    non_cluster_mask = ~cluster_mask
    
    for cell_type, markers in marker_dict.items():
        valid_markers = [m for m in markers if m in adata.var_names]
        
        if len(valid_markers) >= 5:
            mean_in = adata[cluster_mask, valid_markers].X.mean()
            mean_out = adata[non_cluster_mask, valid_markers].X.mean()
            fold_change_scores.loc[cluster, cell_type] = mean_in / mean_out if mean_out > 0 else 0

# Assign based on threshold
FC_THRESHOLD = 1.01
cell_type_assignments = {}

for cluster in clusters:
    top_fc = fold_change_scores.loc[cluster].max()
    if top_fc > FC_THRESHOLD:
        cell_type_assignments[cluster] = fold_change_scores.loc[cluster].idxmax()
    else:
        cell_type_assignments[cluster] = "Unknown"

adata.obs['cell_type'] = adata.obs['cluster'].map(cell_type_assignments)
```

---

## AtlasXOmics-Specific Notes

**Data type:** Gene activity scores from ATAC-seq (not RNA expression)
- Scores typically 5-10% lower than RNA-seq
- Best markers: TFs, surface receptors, metabolic enzymes with regulated expression
- Avoid: secreted proteins, housekeeping genes, post-transcriptional regulators

**Expected structure:**
- `adata.X`: gene activity matrix
- `adata.obs['cluster']`: cluster assignments
- `adata.obsm['spatial']`: spatial coordinates

---

## Common Issues

**All FC values near 1.0 (0.95-1.05):**
- Check marker quality
- May indicate clusters don't represent distinct cell types
- Consider if tissue composition is uniform

**Unknown clusters:**
- First check quality metrics (fragment count, TSS enrichment)
- Extract top 20 genes to identify unexpected cell types
- May represent quiescent/G0 phase if quality is good

**FC > threshold but low separation:**
- Consider "Mixed (A, B)" labels
- Flag for manual review

---

## Validation Requirements

✓ FC > threshold for all assigned cell types  
✓ Unknown clusters explained (quality, quiescent state, or rare type)  
✓ Spatial distribution matches expected anatomy  
✓ Marker expression validated in `w_h5` viewer  
✓ **ALL evaluation metrics computed** (see `evals.md`)  

---

## Key Principles

1. **Never omit expected cell types** — supplement with literature markers if needed
2. **FC > 1 required for assignment** — enrichment, not just presence
3. **Mark weak enrichment as Unknown** — don't force assignments
4. **Document threshold and rationale** — for reproducibility
5. **⚠️ MUST compute evaluation metrics** — read and follow `evals.md`
