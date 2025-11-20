## 10X Xenium Analysis Guideline

This is the **authoritative step-by-step pipeline** for 10X Xenium experiment. Follow steps in order.

---

## Pipeline Overview

1. **Data Loading** — Print analysis header and load data using **Scanpy** and visualize with `w_h5`.
2. **Preprocessing** — Run QC → normalization → PCA → Harmony batch correction → UMAP → Leiden **only if the user confirms**.
3. **Differential Gene Expression (DGE)** — identify marker genes per cluster using rank-based DGE tests (use t-tests by default, but allow the user to select different methods). Report top marker genes for each cluster and make dot plots with scanpy.
4. **Cell Type Annotation** — Aassign biological meaning to clusters using marker-gene dictionaries and provide marker rationales.
5. **Cell Segmentation** — Ask the user whether to resegment. If confirmed, launch the workflow form by running `w_workflow(wf_name="wf.__init__.xenium_cell_segmentation_workflow", ...)`.
6. **Domain Detection** — Convert the final `AnnData` object to `.h5ad` and launch `w_workflow(wf_name="wf.__init__.", ...)`.

---

## General Workflow Rules

- **Always render a form** for any workflow with `lplots.widgets`.
  - Pre-populate sensible defaults.
  - **End with `w.value`** to execute.
- **Summarize results** (and progress/status) using `w_text_output(content=..., appearance={"message_box": “success” | “info”})`.
- **Render figures** using `w_plot` (and table widgets where appropriate). **Do not reuse the variable name fig. Use descriptive names: fig_qc, fig_umap, etc.**
  - Non-Scanpy figures: use **Plotly**.
- Do **not** ask users for file paths or node IDs; auto-discover from the attached folder.
- Keep objects and outputs reproducible: never overwrite without user consent; create new keys/versions.
- **Always include `w.value` at the end** when using a widget.
- Do **not** use `lplots.widgets.number`

---

## Step 1 — Data Loading

**Assumption:** Each sample already has a single, ready-to-use `.h5ad` containing counts, metadata, and spatial coordinates (either adata.obs[Spatial] or adata.obs[X_spatial]).

- Load the `.h5ad` and render with `w_h5(adata)`. If a `spatial_dir` subfolder exists, use `w_h5(adata, spatial_dir)` instead.
- If embeddings/clusters already exist, do not run preprocessing yet. Ask for confirmation in **Step 3** before recomputing.

---

## Step 2 — Preprocessing (run only if the user confirms)

**Goal:** Create a clean, integrated embedding and clusters ready for DGE and annotation.

- **Quality control**
  - Compute standard QC metrics (genes per cell, counts per cell, mitochondrial fraction if available).
    ```python
    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
    ```
  - Show histograms using plotly; expose filter thresholds via widgets `w_text_input`.
  - Confirm before applying filters; keep an unfiltered copy.
  - **ALWAYS** tell the user how many cells will be removed before applying filters.
  - If all cells are removed, ask user to re-enter parameters.

- **Normalization**
  - Library-size normalization, then log transform.

- **Dimensionality reduction**
  - PCA on all genes; store in `adata.obsm["X_pca"]`.
  ```python
    sc.tl.pca(adata, n_comps=n_pcs)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=n_pcs, show=False)

    # Capture the figure and render with the plot widget
    fig_pca = plt.gcf()
    w_plot(label="PCA Explained Variance Ratio", source=fig_pca)

    plt.close()
  ```

- **Neighborhood graph, UMAP, clustering**
  - Build neighbors on the integrated space.
  - Compute UMAP; store in `adata.obsm["X_umap"]`.
  - Run Leiden clustering; store categorical labels in adata.obs["leiden"]; **update the original w_h5 viewer in place—do not create a new viewer.**
  - Render Leiden clusters on the spatial embedding and UMAP using Plotly subplots.

- **Non-destructive policy**
  - Do not overwrite existing embeddings or clusters unless the user opts in. Write to new keys if needed.

### Parameters (render as widgets with sensible defaults)
- **QC thresholds:** `min_genes` (default: 5), `min_counts` (default: 20), `max_mt_frac` (default: 0.2)
- **PCA:** `n_pcs` (default: 50)
- **Neighbors:** `n_neighbors` (default: 15), `metric` (default: `"euclidean"`)
- **Leiden:** `resolution` (default: 0.5), `random_state` (default: 0)

> Generate widgets for the parameters above; pre-fill with defaults. Require explicit user confirmation before running.

## Step 3 — Differential Gene Expression (DGE)

- **Default test (fastest):** `t-test_overestim_var ` with Benjamini–Hochberg FDR control.
- **Alternative (more robust to outliers/zero inflation):**
- **Default custering column:** `K=5`.
> **Heuristic:** keep `t-test_overestim_var ` as default (“Quick mode”). Offer a toggle to “Robust mode” → `wilcoxon` when users prioritize robustness over speed.

### Parameters (render as widgets with defaults)
- **Method:** `t-test_overestim_var ` *(default)* | `wilcoxon`
- **Top genes per cluster:** `top_n` *(default:5)*

```python
valid_cells = ~adata.obs[cluster_col].isna()
n_removed = (~valid_cells).sum()
adata_filtered = adata[valid_cells].copy()
adata_filtered.obs[col] = adata_filtered.obs[col].apply(lambda x: str(x) if not isinstance(x, str) else x)
sc.tl.rank_genes_groups(
    adata_filtered,
    groupby=cluster_col,
    method=dge_method,
    use_raw=False,
    corr_method='benjamini-hochberg'
)
```

### Reporting
- Per cluster: gene, log fold change, p-value, **FDR**.
- Store results under `adata.uns["dge"]`.
- Report top marker genes for each cluster and Make dot plots with scanpy.
- Select the top four biologically meaningful marker genes and color the spatial embedding by their expression in four subplots. Explain why they are biologically meaningful.

```python
sc.pl.dotplot(
    adata,
    var_names=all_marker_genes,
    groupby=cluster_col_name,
    standard_scale='var',
    show=False
)

dotplot_fig = plt.gcf()
w_plot(label="Marker Gene Expression Dot Plot", source=dotplot_fig)
plt.close()
```
---

## Step 4 — Cell Type Annotation

- Assign biological meaning to clusters using marker-gene dictionaries to identify cell types.
- Produce:
  - `adata.obs["cell_type"]` with human-readable labels.
  - A brief **marker rationale** per label (the key genes actually used).
  - A confidence score (e.g., based on marker enrichment consistency).
- Summarize assignments and **update the original `w_h5 viewer` with new cell types in place—do not create a new viewer.**
- Visualize predicted cell-type percentages with a Plotly bar chart.

---

## Step 5 — Cell Segmentation (Optional)

- Ask: “Re-segment cells?”
- If **yes**:
  - **Use the full-resolution morphology image only**: `morphology.ome.tif` (never any file containing “mip”).
  - **Render a form** and launch
    `w_workflow(wf_name="wf.__init__.xenium_cell_segmentation_workflow", …).value`
  - **Expose parameters as widgets** (at minimum: `output_dir`, `level`) with sensible defaults.
- If **no**, proceed without changes.

```python
from lplots.widgets.workflow import w_workflow

params = {
    "transcript": LatchFile("latch:///transcripts.parquet"),
    "tiff": LatchFile("latch:///morphology.ome.tif"),
    "run_name": "my_run",
    "output_dir": LatchDir("latch:///cell_resegmentation_output"),
    "level": "5",
}

w = w_workflow(
    wf_name="wf.__init__.xenium_cell_segmentation_workflow",
    version=None,
    automatic=True,
    label="Cell Segmentation Workflow",
    params=params
)
execution = w.value

if execution is not None:
  res = await execution.wait()

  if res is not None and res.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
      # inspect workflow outputs for downstream analysis
      workflow_outputs = list(res.output.values())
```

---

## Step 6 — Domain Detection

- Rename adata.obsm["Spatial"] to adata.obsm["X_spatial"], then save the processed object to `.h5ad`.
    ```python
    if 'Spatial' in adata.obsm.keys():
        adata.obsm['X_spatial'] = adata.obsm['Spatial']

    # Save the processed H5AD file
    output_h5ad_path = "/tmp/xenium_processed_for_domain_detection.h5ad"
    adata.write_h5ad(output_h5ad_path)
    output_lpath = LPath("latch:///xenium/analysis_output/xenium_processed_for_domain_detection.h5ad")
    output_lpath.upload_from(output_h5ad_path)
    local_path = Path(output_h5ad_path)
    output_lpath.upload_from(local_path)
    ```
- Use this file as input for the domain-detection workflow.

```python
params = {
    "input_file": LatchFile("latch:///xenium_processed_for_domain_detection.h5ad"),
    "run_name": "my_run",
    "output_dir": LatchDir("latch:///Domain_detection_output"),
    "lambda_list": "0.5",
}

w = w_workflow(
    wf_name="wf.__init__.domain_detection_wf",
    version=None,
    automatic=True,
    label="Domain Detection Workflow",
    params=params,
)
execution = w.value

if execution is not None:
  res = await execution.wait()

  if res is not None and res.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
      # inspect workflow outputs for downstream analysis
      workflow_outputs = list(res.output.values())
```

