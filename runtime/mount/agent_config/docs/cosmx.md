## CosMX Spatial RNA-seq Analysis Workflow

This is the **authoritative step-by-step pipeline** for **CosMX** spatial transcriptomics experiments. Follow steps in order.

1. **Data Loading** — load data using **Scanpy** and display with `w_h5`.
2. **Preprocessing** — normalize and log-transform expression data.
3. **Dimensionality Reduction & Clustering** — compute PCA, UMAP, and Leiden clusters using **Scanpy**.
4. **Differential Gene Expression & GSEA** — identify marker genes per cluster and perform gene set enrichment analysis.
5. **Cell Type Annotation** — Use CellGuide marker database (see {marker_gene_annotation_docs})
6. **Manual Validation** (optional) — validate or refine automated annotations when needed.

The sections below define detailed guidelines for each step.

---

### **Data Loading**

CosMX data can be provided in two formats:

1. **Spatial directory** (recommended) — contains `transcripts.h5ad` plus tiled images for high-resolution visualization
2. **Direct H5AD file** — a standalone AnnData file with gene expression and spatial coordinates

#### Option 1: Load from Spatial Directory

- The spatial directory contains:
  - `transcripts.h5ad` — the AnnData file with gene expression and spatial coordinates
  - Tiled images for visualization in the `w_h5` viewer
- Use `LPath` or direct file path to load the spatial directory. **Always** prefer **full, human-readable latch:// paths** (not node IDs).
- The H5AD file is located at `spatial_dir/transcripts.h5ad`.

#### Option 2: Load H5AD File Directly

- Locate and load the H5AD file directly using Scanpy.
- The file contains spatial coordinates in `adata.obsm['spatial']`.

#### Example Implementation

```python
import scanpy as sc
from latch.ldata.path import LPath
from pathlib import Path

# Option 1A: Load from spatial directory (remote Latch data)
spatial_dir = LPath("latch:///path/to/cosmx_spatial_dir")
h5ad_path = spatial_dir / "transcripts.h5ad"
adata = sc.read_h5ad(h5ad_path.path)

# Option 1B: Load from spatial directory (local path)
spatial_dir_local = Path("/path/to/cosmx_spatial_dir")
adata = sc.read_h5ad(spatial_dir_local / "transcripts.h5ad")

# Option 2A: Load H5AD file directly (remote)
adata = sc.read_h5ad("latch:///path/to/cosmx_data.h5ad")

# Option 2B: Load H5AD file directly (local)
adata = sc.read_h5ad("/path/to/cosmx_data.h5ad")

# Display with w_h5
from lplots.widgets.h5 import w_h5
viewer = w_h5(ann_data=adata)
```

---

### **Preprocessing**

CosMX data typically contains a small panel of ~900 genes (targeted panel). Therefore, **do not perform feature selection**.

#### Steps:

1. **Normalization** — normalize each cell to the same total count, then log-transform.
   - Use `sc.pp.normalize_total()` followed by `sc.pp.log1p()`.
   - Expose normalization parameters via widgets (target sum, default: 10,000).
2. **Save raw counts** — store raw counts in `adata.raw` before normalization for downstream analysis.

#### Example Implementation

```python
import scanpy as sc
from lplots.widgets.number import w_number_input

# Widget for target sum
target_sum = w_number_input(
    label="Target sum for normalization",
    default=10000,
    min_value=1000,
    max_value=100000
)

# Store raw counts
adata.raw = adata

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=target_sum.value)
sc.pp.log1p(adata)
```

---

### **Dimensionality Reduction & Clustering**

Perform standard dimensionality reduction and clustering using Scanpy.

#### Steps:

1. **Check for existing embeddings** — check if PCA and UMAP already exist in `adata.obsm`. If so, give the user the option to recompute.

2. **PCA** — compute principal components.

   - Expose parameter: number of components (default: 30).
   - Use `sc.tl.pca()`.

3. **Compute neighbors** — build a k-nearest neighbor graph.

   - Expose parameter: number of neighbors (default: 15).
   - Use `sc.pp.neighbors()`.

4. **UMAP** — compute UMAP embedding for visualization.

   - Use `sc.tl.umap()`.

5. **Leiden clustering** — apply Leiden algorithm for community detection.

   - Expose parameter: resolution (default: 0.5).
   - Use `sc.tl.leiden()`.

6. **Visualization** — always visualize clusters on both UMAP and spatial embeddings using `w_h5`.

#### Example Implementation

```python
import scanpy as sc
from lplots.widgets.number import w_number_input
from lplots.widgets.h5 import w_h5

# Check if embeddings exist
has_pca = 'X_pca' in adata.obsm
has_umap = 'X_umap' in adata.obsm

if has_pca and has_umap:
    from lplots.widgets.select import w_select
    recompute = w_select(
        label="PCA and UMAP already exist. Recompute?",
        options=["No", "Yes"],
        default="No"
    )
    should_compute = recompute.value == "Yes"
else:
    should_compute = True

if should_compute:
    # Widget for PCA components
    n_comps = w_number_input(
        label="Number of PCA components",
        default=30,
        min_value=10,
        max_value=100
    )

    # Widget for neighbors
    n_neighbors = w_number_input(
        label="Number of neighbors",
        default=15,
        min_value=5,
        max_value=100
    )

    # Compute PCA
    sc.tl.pca(adata, n_comps=n_comps.value)

    # Compute neighbors
    sc.pp.neighbors(adata, n_neighbors=n_neighbors.value, n_pcs=n_comps.value)

    # Compute UMAP
    sc.tl.umap(adata)

# Widget for clustering resolution
resolution = w_number_input(
    label="Leiden clustering resolution",
    default=0.5,
    min_value=0.1,
    max_value=2.0,
    step=0.1
)

# Compute neighbors if they don't exist
if 'neighbors' not in adata.uns:
    sc.pp.neighbors(adata, n_neighbors=15)

# Leiden clustering
sc.tl.leiden(adata, resolution=resolution.value)

# Visualize clusters
viewer = w_h5(ann_data=adata)
```

---

### **Differential Gene Expression & GSEA**

Identify marker genes for each cluster and perform gene set enrichment analysis (GSEA) to understand biological pathways.

#### Steps:

1. **Differential Gene Expression** — use `sc.tl.rank_genes_groups()` to identify marker genes per cluster.

   - Default method: `wilcoxon` (robust for scRNA-seq data).
   - Compare each cluster against all others (cluster vs. rest).

2. **Visualize top markers** — display top marker genes using `w_table()` and generate dot plots.

3. **Automated Cell Type Annotation** — **immediately after** computing marker genes, automatically run the marker gene-based annotation workflow to suggest cell types. See {marker_gene_annotation_docs} for complete details. This provides data-driven cell type suggestions before manual validation.

4. **Gene Set Enrichment Analysis (GSEA)** — run GSEA on ranked gene lists for each cluster.
   - Use `gseapy.prerank()` for preranked GSEA.
   - Gene sets: use Enrichr libraries (e.g., `KEGG_2021_Human`, `Hallmark_2020`, `Reactome_2022`).
   - Ranking metric: use test statistic scores or log fold changes.

#### Example Implementation

```python
import scanpy as sc
import os
import numpy as np
import pandas as pd
import gseapy as gp
from lplots.widgets.table import w_table
from lplots.widgets.select import w_select

# Differential gene expression
sc.tl.rank_genes_groups(
    adata,
    groupby='leiden',
    method='wilcoxon',
    use_raw=True,
    pts=True,
    tie_correct=True
)

# Display top marker genes for each cluster
for cluster in adata.obs['leiden'].cat.categories:
    marker_df = sc.get.rank_genes_groups_df(adata, group=cluster)
    marker_df = marker_df.head(20)
    w_table(source=marker_df, label=f"Top 20 markers for cluster {cluster}")

# Generate dot plot
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)

# GSEA on Leiden clusters
def gsea_on_leiden(
    adata,
    groupby="leiden",
    gene_sets="KEGG_2021_Human",
    rank_by="scores",
    outdir="gsea_results",
    n_permutations=1000,
    method="wilcoxon",
    use_raw=True,
    min_size=10,
    max_size=5000,
    seed=0,
):
    """
    For each Leiden cluster, compute a cluster-vs-rest ranking (via scanpy.tl.rank_genes_groups),
    then run GSEApy.prerank on the ranked list.

    Returns: dict {cluster_label: (res2d pandas DataFrame, GSEApy result object)}
    """
    np.random.seed(seed)
    os.makedirs(outdir, exist_ok=True)

    if use_raw and adata.raw is not None:
        Xref = adata.raw.to_adata()
        # align obs/var with adata
        Xref = Xref[:, Xref.var_names.intersection(adata.var_names)]
        ad = adata.copy()
        ad.X = ad[:, Xref.var_names].X
        ad.var = Xref.var.loc[ad.var_names]
    else:
        ad = adata

    if groupby not in ad.obs:
        raise ValueError(f"'{groupby}' not found in adata.obs")

    # Differential ranking (cluster vs rest) to get a statistic for GSEA ranking
    sc.tl.rank_genes_groups(
        ad,
        groupby=groupby,
        method=method,
        use_raw=False,
        pts=True,
        tie_correct=True
    )

    clusters = ad.obs[groupby].astype(str).unique().tolist()
    results = {}

    for cl in clusters:
        # Extract per-cluster table
        df = sc.get.rank_genes_groups_df(ad, group=cl)
        # Choose ranking metric
        if rank_by not in df.columns:
            raise ValueError(f"rank_by='{rank_by}' not available. Options include: {list(df.columns)}")
        # Build ranked list (gene, score), descending
        rnk = df.loc[:, ["names", rank_by]].dropna()
        # Avoid duplicate genes (keep the best-ranked)
        rnk = rnk.sort_values(rank_by, ascending=False).drop_duplicates("names", keep="first")
        rnk.columns = ["gene", "score"]

        # Run prerank GSEA
        cl_out = os.path.join(outdir, f"{groupby}_{cl}")
        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets=gene_sets,
            outdir=cl_out,
            min_size=min_size,
            max_size=max_size,
            permutation_num=n_permutations,
            seed=seed,
            format="png",
            no_plot=False
        )
        # Save the ranked list used
        rnk.to_csv(os.path.join(cl_out, "ranked_list.tsv"), sep="\t", index=False)
        # Collect table of enriched sets (NES, pval, FDR)
        results[cl] = (pre_res.res2d, pre_res)

    return results

# Run GSEA
gsea_results = gsea_on_leiden(
    adata,
    groupby="leiden",
    gene_sets="KEGG_2021_Human",
    rank_by="scores",
    outdir="gsea_results",
    n_permutations=1000,
    method="wilcoxon",
    use_raw=True,
    seed=0
)

# Display GSEA results for each cluster
for cluster, (res_df, res_obj) in gsea_results.items():
    top_pathways = res_df.head(10)
    w_table(source=top_pathways, label=f"Top 10 enriched pathways for cluster {cluster}")
```

---

### **Cell Type Annotation**

Assign biological cell type identities to clusters based on spatial data.

#### CellGuide Marker Gene Annotation (Any Organism)

Use CellGuide marker gene database annotation for human or mouse data.

See {marker_gene_annotation_docs} for the complete workflow.

**Brief overview:**

1. Load CellGuide marker gene databases (per-gene and per-celltype)
2. Confirm organism and tissue type with the user
3. For each cluster, query the per-gene database with top 10 marker genes
4. Aggregate cell type suggestions across markers using consensus scoring
5. Validate with per-celltype database (canonical markers)
6. Display suggested cell types with confidence scores
7. Visualize on UMAP and spatial embeddings using `w_h5`

**Advantages:**

- Works for any organism (mouse, human, etc.)
- Fast (< 1 minute)
- Discrete, interpretable labels
- Transparent evidence trail

**Limitations:**

- Discrete labels only (no probabilities)
- Limited to database coverage
- May miss novel cell types

---

#### Manual Validation (Optional)

Use this phase to:

- **Validate** automated suggestions when confidence is low
- **Refine** annotations with additional known markers
- **Discover** cell subtypes not well-represented in the database

**Steps:**

1. **Review automated results** — check the suggested cell types and confidence scores from Phase 1.

2. **Identify clusters needing validation** — focus on:

   - Low-confidence predictions
   - "Unknown" annotations
   - Unexpected cell types

3. **Define validation markers** — for clusters requiring manual review, specify additional marker genes:

   - Use domain knowledge for the specific tissue
   - Include canonical markers from literature
   - Consider both positive and negative markers

4. **Compute gene set scores** — use `sc.tl.score_genes()` to compute marker enrichment scores.

5. **Compare with automated suggestions** — visualize both automated predictions and validation scores using `w_h5` and Plotly.

#### Example Implementation: Phase 2 (Manual Validation)

**Note:** Phase 1 (automated annotation) runs automatically after differential expression. The code below is for **optional manual validation** only.

```python
import scanpy as sc
from lplots.widgets.text import w_text_input
from lplots.widgets.select import w_select
from lplots.widgets.h5 import w_h5
from lplots.widgets.table import w_table
import pandas as pd

# Review automated results (from Phase 1)
# Check which clusters have low confidence or unexpected annotations
low_confidence_clusters = adata.obs[
    adata.obs['annotation_confidence'] < 10
]['leiden'].unique()

if len(low_confidence_clusters) > 0:
    # Manual validation for low-confidence clusters

    # Define validation marker genes (example for brain tissue)
    validation_markers = {
        "Excitatory Neurons": ["SLC17A7", "NEUROD6", "SATB2"],
        "Inhibitory Neurons": ["GAD1", "GAD2", "SLC32A1"],
        "Astrocytes": ["GFAP", "AQP4", "SLC1A2"],
        "Oligodendrocytes": ["MBP", "MOG", "PLP1"],
        "Microglia": ["C1QA", "C1QB", "CSF1R"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5"]
    }

    # Compute validation scores
    for cell_type, markers in validation_markers.items():
        markers_present = [m for m in markers if m in adata.var_names]
        if markers_present:
            score_name = f"{cell_type}_validation_score"
            sc.tl.score_genes(adata, markers_present, score_name=score_name)

    # Compare automated vs validation results
    comparison_data = []
    for cluster in low_confidence_clusters:
        cluster_mask = adata.obs['leiden'] == cluster
        predicted = adata.obs.loc[cluster_mask, 'predicted_cell_type'].iloc[0]
        auto_conf = adata.obs.loc[cluster_mask, 'annotation_confidence'].iloc[0]

        # Get validation scores
        validation_cols = [col for col in adata.obs.columns if col.endswith("_validation_score")]
        if validation_cols:
            mean_scores = adata.obs.loc[cluster_mask, validation_cols].mean()
            top_validation = mean_scores.idxmax().replace("_validation_score", "")
            top_val_score = mean_scores.max()
        else:
            top_validation = "N/A"
            top_val_score = 0

        comparison_data.append({
            "Cluster": cluster,
            "Automated Prediction": predicted,
            "Auto Confidence": auto_conf,
            "Manual Validation": top_validation,
            "Validation Score": top_val_score
        })

    comparison_df = pd.DataFrame(comparison_data)
    w_table(source=comparison_df, label="Automated vs Manual Validation Comparison")

    # Visualize with w_h5 to inspect both predictions
    viewer = w_h5(ann_data=adata)
else:
    # All clusters have high confidence
    print("All clusters have high-confidence automated annotations. Manual validation not needed.")
```

---

## ✅ **General Rules**

- **DO NOT** delete cells unless explicitly prompted by the user.
- **DO NOT** create duplicated cells.
- Always use markdown to summarize text output after each analysis step. **DO NOT** use `print()`.
- Always produce **visual outputs** at every stage. Use `w_plot()` and `w_table()` to display plots and tables on the UI.
- Always use `w_h5` for spatial embedding and UMAP. For all other plots, use Plotly instead. Do not use both `w_h5` and Plotly for the same visualization to avoid redundancy.
- Always pass a DataFrame object directly to the `source` argument of `w_table`. Do not pass a method call or expression (e.g., `df.head()`, `df.round(3)`, `df.sort_values(...)`) directly into `w_table`.
- Apply sensible defaults for all parameters.
- Add well-formatted markdown summaries.
- Prioritize **transparency**, **reproducibility**, and **interactivity** throughout the analysis.

---

## Data Assumptions

CosMX datasets can be provided in two formats:

### Format 1: Spatial Directory (Recommended)

A directory containing:

- `transcripts.h5ad` — the main AnnData file with gene expression and spatial coordinates
- Tiled images for high-resolution visualization in the `w_h5` viewer

#### Spatial Directory Structure

```
spatial_dir/
├── transcripts.h5ad
└── (other spatial data files like images and coordinates, etc)
```

### Format 2: Direct H5AD File

A standalone H5AD file containing gene expression and spatial coordinates.

### H5AD File Contents

Both formats contain the same H5AD structure:

- **Gene expression matrix**: raw or normalized counts for ~900 targeted genes.
- **Spatial coordinates**: stored in `adata.obsm['spatial']`.
- **Cell-level metadata**: stored in `adata.obs` (may include sample, condition, FOV, cell_id, etc.).

Standard AnnData structure:

```
AnnData object with n_obs × n_vars = 50000 × 900
    obs: 'sample', 'fov', 'cell_id'
    var: 'gene_name'
    obsm: 'spatial'
```

### Loading Guidelines

- **Always** assume spatial coordinates are in `adata.obsm['spatial']`.
- If provided a spatial directory, load data from `spatial_dir/transcripts.h5ad`.
- If provided a direct H5AD file, load it directly.
- Do not assume the presence of pre-computed PCA, UMAP, or clusters unless verified.
