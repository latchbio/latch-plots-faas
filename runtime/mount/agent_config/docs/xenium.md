## Analysis Guideline

This is the **authoritative step-by-step pipeline** for 10X Xenium experiment. Follow steps in order.

1. **Data Loading** — Load data using **Scanpy** and visualize with `w_h5`.
2. **Cell Type Annotation** — This step uses **precomputed differential expression results** from the `analysis/diffexp/` directory.
The agent must **not** perform new differential expression analysis unless specified by user.
3. **Cell Segmentation** — Ask the user whether to resegment. If confirmed, run `w_workflow(wf_name="wf.__init__.xenium_cell_segmentation_workflow", ...)`.
4. **Domain Detection** — Convert the final `AnnData` object to `.h5ad` and launch `w_workflow(wf_name="wf.__init__.domain_detection_wf", ...)`.


The section below defines detailed guidelines for each of the above steps.

### **General Workflow Rules**
When the step requires launching a workflow:
- **Always render a form** with `lplots.widgets` for all workflow inputs. **Do not** ask users questions in text. Pre-populate each widget with sensible `default` whenever possible.
- **Always include `w.value` at the end**
- **Always summarize results using** `w_text_output`
- **Always render figures with** `w_plot`

### **Progress Messages Policy**
When describing the current step or action (e.g., “Loading precomputed differential expression results…” or “Using known marker genes to identify major brain cell types…”), always use the info message box appearance.

```python
from lplots.widgets.text import w_text_output
w_text_output(
    content="Loading precomputed differential expression results...",
    appearance={"message_box": "info"}
)
```
---

## Step 1 — Data Loading

#### **Case A — Single Xenium Output Folder**
The user will attach a **single folder** containing all Xenium outputs (e.g. `cell_feature_matrix.h5`, `cells.csv.gz`, `analysis/`, etc.).
All files and subfolders (e.g. `umap/`, `pca/`, `clustering/`, `diffexp/`) must be located **automatically** by the agent — do **not** ask for file paths or node IDs.

Use the attached folder as the root directory for all subsequent file loading.

#### **Case B — Two Xenium Output Folders**
The user will attach **two folders**, each representing a separate Xenium output (e.g. `cell_feature_matrix.h5`, `cells.csv.gz`, `analysis/`, etc.).
All files and subfolders (e.g. `umap/`, `pca/`, `clustering/`, `diffexp/`) must be located **automatically** by the agent — do **not** ask for file paths or node IDs.

For each folder:
- Perform the same loading steps (**1.1–1.5**) as described below.
- Add metadata such as `sample_id` and any relevant identifiers to `adata.obs` for that sample.
- Store each sample’s data in its own `AnnData` object.
- **After loading each sample**, compute and visualize **quality-control metrics** in a single figure:
  - The figure must contain **two subplots** side by side rendered with `w_plot`:
    1. Histogram of **total transcripts per cell** (`total_counts`)
    2. Histogram of **number of detected genes per cell** (non-zero features in `adata.X`)
  - Report the **number of cells**, **number of genes**, and the **median** and **interquartile range (IQR)** of both metrics using `w_text_output`.

After both samples are loaded:
1. **Do not merge** the AnnData objects — keep each sample separate.
2. **Render a selection widget** (e.g. `w_select`) that allows the user to choose which AnnData object to visualize. The widget should clearly display sample identifiers (e.g. folder names or sample IDs) and update the visualization dynamically based on the selected sample.


- The loaded Xenium data (expression matrix, spatial coordinates, embeddings, and clustering) is **already preprocessed** by 10x Genomics pipelines.
- The agent must **not** perform additional single-cell preprocessing steps such as:
  - Quality control or filtering
  - Normalization
  - Feature selection
  - Re-running PCA, UMAP, or clustering
- These steps should **only** be re-executed if the user **explicitly requests** reprocessing (for example, to change parameters or regenerate embeddings).
- By default, assume the loaded AnnData object is **analysis-ready**.

### 1.1 Locate Files and Extract Analysis Archive

Search recursively within the attached folder for:

| File / Folder | Purpose |
|----------------|----------|
| `cell_feature_matrix.h5` | Feature matrix input |
| `cells.csv.gz` | Cell metadata and spatial coordinates |
| `analysis/umap/...` | UMAP embeddings |
| `analysis/pca/...` | PCA embeddings |
| `analysis/tsne/...` | t-SNE embeddings |
| `analysis/clustering/...` | Cluster assignments |
| `analysis/diffexp/...` | Differential-expression results |

If the file `analysis.tar.gz` is found, **always automatically extract it** into a new subdirectory named `analysis/` before proceeding.
Do not ask the user to confirm extraction.

```python
import tarfile, os

analysis_tar = os.path.join("<attached_folder>", "analysis.tar.gz")
if os.path.exists(analysis_tar):
    with tarfile.open(analysis_tar, "r:gz") as tar:
        tar.extractall(path="<attached_folder>")
```

### 1.2 Load Feature Matrix
```python
adata = read_10x_h5(LPath("<attached_folder>/cell_feature_matrix.h5"))
```

### 1.3 Load Spatial Coordinates
```python
import pandas as pd
df = pd.read_csv("<attached_folder>/cells.csv.gz")

adata.obsm["spatial"] = df[["x_centroid", "y_centroid"]].to_numpy()
adata.obs["transcript_counts"] = df["transcript_counts"]
adata.obs["total_counts"] = df["total_counts"]
```

### 1.4 Load Analysis Outputs

Extract data from the `analysis/` subfolder and attach them to the existing `AnnData` object:

| File | Target Field | Description |
|------|---------------|-------------|
| `analysis/umap/gene_expression_2_components/projection.csv` | `adata.obsm["X_umap"]` | UMAP coordinates |
| `analysis/pca/gene_expression_10_components/projection.csv` | `adata.obsm["X_pca"]` | PCA coordinates |
| `analysis/tsne/gene_expression_2_components/projection.csv` | `adata.obsm["X_tsne"]` | t-SNE coordinates |
| `analysis/clustering/gene_expression_kmeans_x_clusters/clusters.csv` | `adata.obs["cluster_x"]` | Cluster assignments (categorical), load as string |
| `analysis/diffexp/gene_expression_kmeans_x_clusters/differential_expression.csv` | `adata.uns["diffexp_genes"]` | Differentially expressed genes (from column `Feature Name`) |

**Never** assign UMAP, PCA, or t-SNE arrays directly from CSV files. **Always** use `set_index` and `reindex` to match embeddings and clustering data to the main AnnData object.

### 1.5 Display the Final Object
Render the merged object using:
```python
w_h5(adata)
```

---

## Step 2 — Cell Type Annotation

This step uses **precomputed differential expression results** from the `analysis/diffexp/` directory.
The agent must **not** perform new differential-expression analysis unless user specifies.


### Cluster Label Normalization
Before using any `cluster_*` column, standardize labels to **strings** and mark missing as `"unassigned"`.
Remove trailing decimals so `"1.0" → "1"`, `"2.0" → "2"`. Never compare numerically.

```python
adata.obs[cluster_col] = (
    adata.obs[cluster_col]
    .astype(str)
    .str.replace(r"\.0$", "", regex=True)
    .fillna("unassigned")
    .astype("category")
)
```

### 2.1 Ask User for Cluster Selection
- Render a widget that allows the user to select the desired number of clusters (2–10):

```python
from lplots.widgets import w_select

w = w_select(
    label="Select number of clusters",
    options=[str(i) for i in range(2, 11)],
    default="5"
)
n_clusters = int(w.value)

```

### 2.2 Locate and Load Differential Expression Results
**Always use this exact relative path. Do not alter or guess another directory**.
- **Path**: analysis/diffexp/gene_expression_kmeans_<k>_clusters/differential_expression.csv
- **Display** top markers statistics from differential_expression.csv in a table using `w_table`.

```python
diffexp_path = f"analysis/diffexp/gene_expression_kmeans_{n_clusters}_clusters/differential_expression.csv"
df_diffexp = pd.read_csv(f"<attached_folder>/{diffexp_path}")
```

### 2.3 Marker Source & Matching
Marker source priority (no improvisation):
1. User-provided markers (if given), case-insensitive.
2. Top markers from `df_diffexp` (e.g., highest fold-change and lowest adjusted p-value):

```python
if diffexp_path.exists():
    w_text_output(content="**Top 5 markers per cluster (by fold change):**")

    cluster_ids = sorted({
        int(c.split()[1])
        for c in df_diffexp_s1.columns
        if c.startswith("Cluster ") and "Log2 fold change" in c
    })

    for cluster_id in cluster_ids:
        fc_col = f"Cluster {cluster_id} Log2 fold change"
        name_col = "Feature Name"

        # Take top 5 genes with highest log2 fold change for this cluster
        cluster_markers = (
            df_diffexp_s1[[name_col, fc_col]]
            .sort_values(fc_col, ascending=False)
            .head(5)
        )
        markers_list = ", ".join(cluster_markers[name_col].tolist())
        w_text_output(content=f"- **Cluster {cluster_id}**: {markers_list}")
```

### 2.4 Derive and Assign Cell Type Labels
- Use the resolved marker set.
- Compute a **confidence score** based on mean expression of top markers. For each cluster, compute the mean expression of the top marker genes across all its cells. Scale these values to the range [0, 1] using min–max normalization.
- Infer biological cell-type names and add to anndata object.

```python
adata.obs["cell_types"] = inferred_labels
adata.obs["ct_confidence"] = confidence_scores
```

- Use `w_text_output` to summarize:
    - Number of annotated clusters
    - Example top markers per cluster
    - Median confidence score

### 2.5 Visualization
- Show a heatmap of expression for the **final chosen markers** per cluster
- Render the updated AnnData object with:

```python
w_h5(adata)
```

---

## Step 3 — Cell Segmentation

- Ask the user: “Do you want to re-segment the cells?”
- If **yes**, launch the segmentation workflow with pre-populated fields:

```python
from lplots.widgets.workflow import w_workflow

params = {
    "transcript": LatchFile("latch:///transcripts.parquet"),
    "tiff": LatchFile("latch:///morphology.ome.tif"),
    "run_name": "my_run",
    "output_dir": LatchDir("latch:///cell_resegmentation_output"),
}

w = w_workflow(
    wf_name="wf.__init__.xenium_cell_segmentation_workflow",
    version=None,
    label="Launch Cell Segmentation Workflow"
)
execution = w.value
```

---

## Step 4 — Domain Detection

- Save the processed object to `.h5ad`.
- Use this file as input for the domain-detection workflow.

```python
params = {
    "input_file": LatchFile("latch:///.h5ad"),
    "run_name": "my_run",
    "output_dir": LatchDir("latch:///Domain_detection_output"),
    "lambda_list": "0.5",
}

w = w_workflow(
    wf_name="wf.__init__.domain_detection_wf",
    version=None,
    label="Launch Domain Detection Workflow"
)
execution = w.value
```

---

