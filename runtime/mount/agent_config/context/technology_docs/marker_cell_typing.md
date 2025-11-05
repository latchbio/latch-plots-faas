# Marker-based Cell Type Annotation

## Overview
Annotate clusters in single-cell or spatial datasets by identifying **marker genes** and matching them to known cell types using the **CellGuide** database.  
This workflow detects marker genes, finds consensus across samples, and performs database lookups to suggest likely cell identities.

---

## Requirements
- **AnnData (`adata`)** containing:
  - `adata.obs["cluster"]` — cluster IDs (required)  
  - Optional sample/condition column: `"sample"`, `"batch"`, `"condition"`, etc.
- **Gene expression data** in `adata.X` or a layer  
- **CellGuide databases:**
  - `cellguide_marker_gene_database_per_gene.json` → annotation lookups  
  - `cellguide_marker_gene_database_per_celltype.json` → validation/reference  
  - `tissues_per_organism.json` → valid tissue options  
- **Organism / tissue context** (e.g., `"Mus musculus"`, `"brain"`)  

Always inspect `adata.obs.columns` before proceeding — column names vary by dataset.

---

## Database Summary

| Database | Structure | Purpose | Example |
|-----------|------------|----------|----------|
| **Per-gene** | `gene → organism → tissue → [cell types]` | Automated annotation | `marker_db_per_gene["Gad1"]["Mus musculus"]["brain"]` → `["inhibitory neuron", "GABAergic neuron"]` |
| **Per-celltype** | `organism → tissue → cell_type → [genes]` | Validation and lookup | — |
| **Tissue options** | `organism → [tissues]` | Populate dropdowns | — |

All databases support organism/tissue-specific lookups, with `"All Tissues"` as fallback.

---

## Pipeline

### Step 0 — Inspect AnnData
Check available metadata columns to identify cluster and sample identifiers.  
If multiple samples exist → run per-sample analysis.  
Otherwise → treat as single-sample dataset.

---

### Step 1 — Differential Expression (Marker Genes)
Use **Wilcoxon rank-sum test** per cluster:

```python
sc.tl.rank_genes_groups(adata, groupby="cluster", method="wilcoxon", n_genes=300)
df = sc.get.rank_genes_groups_df(adata, group=None)
```

If multiple samples exist:
1. Split by sample column  
2. Run `rank_genes_groups` per subset  
3. Concatenate results with `pd.concat()`  

Output columns: `["names", "group", "pvals_adj", "scores", "sample_id"]`

---

### Step 2 — Consensus Markers
Identify reproducible genes across samples:

- Filter where `pvals_adj < 0.05`
- Count across samples
- Define:
  - **`genes_all`** → present in all samples (100%)  
  - **`genes_70`** → present in ≥70% of samples  

Consensus markers are more reliable for annotation.

---

### Step 3 — CellGuide Lookup
Query the database:

```python
db[gene][organism][tissue or "All Tissues"]
```
Fallback to `"All Tissues"` if tissue-specific data missing.

**Inputs:** `organism`, `tissue`, `db_path`.

---

### Step 4 — Summarize Cluster Annotations

Count marker support per cell type and select top hits:

```python
summary[cluster] = {
  "genes_all": genes_all,
  "genes_70": genes_70,
  "most_common_cell_type": top_ct,
  "cell_type_counts": cell_type_counts,
  "markers_for_most_common_cell_type": top_markers,
}
```

---

### Step 5 — Annotate AnnData

Map clusters → cell types:

```python
adata.obs["cell_type"] = (
    adata.obs["cluster"]
    .astype(str)
    .map(mapping)
    .fillna("Unknown")
)
```

Join multiple cell types with commas; assign `"Unknown"` when no match found.

---

## Core Functions

```python
def find_overlap_genes(cluster, df):
    g = df[df["group"] == cluster]
    g = g[g["pvals_adj"] < 0.05]
    if g.empty: return [], []
    overlap = g.groupby("names")["sample_id"].nunique().reset_index(name="n")
    n_total = g["sample_id"].nunique()
    genes_all = overlap.loc[overlap.n == n_total, "names"]
    genes_70 = overlap.loc[overlap.n >= 0.7 * n_total, "names"]
    return list(genes_all), list(genes_70)

def lookup_cellguide_celltypes(genes, organism, tissue, db_path):
    import json, pathlib
    db = json.load(open(pathlib.Path(db_path)))
    res = {}
    for g in genes:
        org = db.get(g, {}).get(organism)
        if not org: continue
        res[g] = org.get(tissue) or org.get("All Tissues") or []
    return res

def summarize_clusters(adata, df, db_path, organism, tissue):
    from collections import Counter
    summary, celltype_markers = {}, {}
    for c in sorted(adata.obs["cluster"].unique()):
        genes_all, genes_70 = find_overlap_genes(c, df)
        if not genes_all:
            summary[c] = {"genes_all": [], "genes_70": genes_70,
                          "most_common_cell_type": [], "cell_type_counts": {},
                          "markers_for_most_common_cell_type": []}
            continue
        lookup = lookup_cellguide_celltypes(genes_all, organism, tissue, db_path)
        counter = Counter(ct for cts in lookup.values() for ct in cts)
        if not counter:
            summary[c] = {"genes_all": genes_all, "genes_70": genes_70,
                          "most_common_cell_type": [], "cell_type_counts": {},
                          "markers_for_most_common_cell_type": []}
            continue
        max_cts = [ct for ct, n in counter.items() if n == max(counter.values())]
        top_markers = [g for g, cts in lookup.items() if any(ct in max_cts for ct in cts)]
        summary[c] = {"genes_all": genes_all, "genes_70": genes_70,
                      "most_common_cell_type": max_cts,
                      "cell_type_counts": dict(counter),
                      "markers_for_most_common_cell_type": top_markers}
        for ct in max_cts:
            celltype_markers.setdefault(ct, set()).update(top_markers)
    return summary, {ct: sorted(list(m)) for ct, m in celltype_markers.items()}
```

---

## Tie-breaking Strategies
If multiple cell types have equal support:
1. Prefer stronger average logFC markers  
2. Choose the more specific ontology term  
3. Favor tissue-relevant cell types  
4. Use `genes_70` for cross-sample consistency  
5. Defer to manual review if still ambiguous  

---

## Practical Notes
- Best for: well-characterized tissues with known markers  
- Not ideal for: rare or novel cell types  
- If empty results: check organism/tissue match or marker thresholds  
- Output: adds `adata.obs["cell_type"]` for downstream visualization  

---
