# **Approach 2 — Cluster-based Annotation**

This approach infers cell types by identifying **marker genes per cluster**, comparing them to known CellGuide markers, and assigning the most likely identity.

## Prerequisites

Your `AnnData` must include:

- `adata.obs["cluster"]` — cluster assignments (required)
- Optionally, one of: `"sample"`, `"condition"`, or `"batch"` — for per-sample analysis
- Gene expression data in `adata.X` or `adata.layers`

Also required:
- **CellGuide database** (`cellguide_marker_gene_database_per_gene.json`)
- **Organism** and **tissue** context (e.g., `"Mus musculus"`, `"brain"`)

> Tip: Always inspect `adata.obs.columns` before making assumptions about metadata column names.

## Step 0 — Inspect AnnData Structure

Identify available metadata columns and determine if the dataset includes multiple samples.

```python
print(adata.obs.columns)
```

If you find a sample column (e.g., `"sample"`, `"condition"`, or `"batch"`) with multiple unique values, perform **per-sample** analysis.  
Otherwise, treat it as a single-sample dataset.

## Step 1 — Differential Expression Analysis

Identify marker genes per cluster using a **Wilcoxon rank-sum test**.  
If multiple samples exist, perform the test within each sample and combine results.

### Single-sample dataset
```python
import scanpy as sc, pandas as pd

PVAL_MAX, LFC_MIN, TOP_N = 0.05, 0.5, 300

sc.tl.rank_genes_groups(adata, groupby="cluster", method="wilcoxon", n_genes=None)
df = sc.get.rank_genes_groups_df(adata, group=None)

df = df[(df.pvals_adj < PVAL_MAX) & (df.logfoldchanges > LFC_MIN)]
df = df.sort_values(["group", "logfoldchanges", "scores"], ascending=[True, False, False])
ranked_genes_df = df.groupby("group").head(TOP_N).reset_index(drop=True)
```

### Multi-sample dataset
```python
all_results = []
for s in adata.obs["sample"].unique():
    sub = adata[adata.obs["sample"] == s].copy()
    sc.tl.rank_genes_groups(sub, groupby="cluster", method="wilcoxon", n_genes=None)
    df = sc.get_rank_genes_groups_df(sub, group=None)
    df["sample"] = s
    df = df[(df.pvals_adj < 0.05) & (df.logfoldchanges > 0.5)]
    df = df.sort_values(["group", "logfoldchanges", "scores"], ascending=[True, False, False])
    all_results.append(df.groupby("group").head(300))
ranked_genes_df = pd.concat(all_results)
```

## Step 2 — Identify Consensus Marker Genes

Determine genes consistently upregulated across samples.

```python
def find_overlap_genes(cluster, ranked_genes_df):
    g = ranked_genes_df[ranked_genes_df["group"] == cluster]
    g = g[g["pvals_adj"] < 0.05]
    if g.empty: return [], []
    overlap = g.groupby("names")["sample"].nunique().reset_index(name="n_conditions")
    n_total = g["sample"].nunique()
    genes_all = overlap.loc[overlap.n_conditions == n_total, "names"]
    genes_70 = overlap.loc[overlap.n_conditions >= 0.7 * n_total, "names"]
    return list(genes_all), list(genes_70)
```

- `genes_all`: markers consistent in **all** conditions (100%)  
- `genes_70`: markers consistent in **≥70%** of conditions

## Step 3 — Lookup Cell Types in CellGuide

Map consensus markers to known cell types.

```python
import json
from collections import Counter
from pathlib import Path

def lookup_cellguide_celltypes(genes, organism, tissue, db_path):
    with open(db_path) as f:
        db = json.load(f)
    results = {}
    for g in genes:
        org_entry = db.get(g, {}).get(organism, {})
        cts = org_entry.get(tissue) or org_entry.get("All Tissues")
        if cts:
            results[g] = cts
    return results
```

## Step 4 — Summarize Cluster Annotations

Combine marker overlap and CellGuide results into a summary dictionary.

```python 
from collections import Counter
import json

def summarize_clusters(adata, ranked_df, db_path, organism="Mus musculus", tissue="brain", min_markers=3):
    with open(db_path) as f:
        db = json.load(f)

    def lookup_cellguide_celltypes(genes):
        results = {}
        for g in genes:
            org_entry = db.get(g, {}).get(organism, {})
            cts = org_entry.get(tissue) or org_entry.get("All Tissues")
            if cts:
                results[g] = cts
        return results

    summary = {}
    for c in sorted(adata.obs["cluster"].unique()):
        genes_all, genes_70 = find_overlap_genes(c, ranked_df)
        lookup = lookup_cellguide_celltypes(genes_all)

        if not lookup:
            summary[c] = {
                "genes_all": genes_all,
                "genes_70": genes_70,
                "most_common_cell_type": [],
                "cell_type_counts": {},
                "markers_for_most_common_cell_type": [],
            }
            continue

        # Count supporting markers per cell type
        counter = Counter(ct for v in lookup.values() for ct in v)

        # Enforce threshold: only consider cell types with ≥ min_markers
        eligible = {ct: n for ct, n in counter.items() if n >= min_markers}
        if not eligible:
            summary[c] = {
                "genes_all": genes_all,
                "genes_70": genes_70,
                "most_common_cell_type": [],
                "cell_type_counts": dict(counter),  # keep raw counts for debugging
                "markers_for_most_common_cell_type": [],
            }
            continue

        max_support = max(eligible.values())
        top_cts = [ct for ct, n in eligible.items() if n == max_support]

        # Markers that support at least one top cell type
        top_markers = [g for g, v in lookup.items() if any(ct in top_cts for ct in v)]

        summary[c] = {
            "genes_all": genes_all,
            "genes_70": genes_70,
            "most_common_cell_type": top_cts,
            "cell_type_counts": dict(counter),
            "markers_for_most_common_cell_type": top_markers,
        }
    return summary
```

**Outputs:**
- `most_common_cell_type`: top predicted identities per cluster  
- `cell_type_counts`: gene support for each type  
- `markers_for_most_common_cell_type`: genes supporting top label(s)

## Step 5 — Annotate `adata`

Add annotations directly to `adata.obs` for downstream plotting.

```python
def annotate_adata(adata, cluster_summary, cluster_col="cluster", cell_col="cell_type"):
    mapping = {
        str(c): ", ".join(v.get("most_common_cell_type", []) or ["Unknown"])
        for c, v in cluster_summary.items()
    }
    adata.obs[cell_col] = adata.obs[cluster_col].astype(str).map(mapping).fillna("Unknown")
    print(f"Annotated {adata.obs[cell_col].nunique()} unique cell types")
    return adata
```

---

## Best Practices and Troubleshooting

| Issue | Likely Cause | Recommendation |
|------|--------------|----------------|
| No CellGuide matches | Tissue or organism mismatch | Verify correct organism/tissue keys |
| Many “Unknown” clusters | Weak markers | Relax p-value or fold-change thresholds |
| Multiple top cell types | Shared marker profiles | Keep both; could represent mixed lineage |
| Missing cluster column | Dataset not clustered | Run UMAP + Leiden/Louvain before annotation |

---

### Next Steps
- Refer to <atlasx_analysis_overview> to see the guide for other types of analyses. 
