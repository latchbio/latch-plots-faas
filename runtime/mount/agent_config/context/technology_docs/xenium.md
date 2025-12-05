## 10X Xenium Analysis Guideline

---

## Pipeline Overview


0. **Data Preparation*** - Convert raw Xenium output to an h5 viewer–compatible format by launching
`w_workflow(wf_name="wf.__init__.xenium_preprocess_workflow", ...)`. **only when there is no `.h5ad` file in the attached data**.
  - If an `.h5ad` is already present in the attached folder, **skip this step** and use that file directly for downstream analysis.
  - If **no `.h5ad` file is present**, **automatically launch workflow without asking the user for confirmation**.
  - **Do NOT treat the following as an existing `.h5ad`**:
    - `analysis.zarr.zip`
    - `analysis.zarr`
    - `analysis.tar.gz`
    - `cell_feature_matrix.h5`
    - Any other generic `.h5` file
    - Any `.zarr` or `.zarr.zip` file
1. **Data Loading** — Print analysis header and load data using **Scanpy** and visualize with `w_h5`.
2. **Preprocessing** — Run QC → normalization → PCA → Harmony batch correction → UMAP → Leiden **only if the user confirms**.
3. **Differential Gene Expression (DGE)** — identify marker genes per cluster using rank-based DGE tests (use t-tests by default, but allow the user to select different methods). Report top marker genes for each cluster and make dot plots with scanpy.
4. **Cell Type Annotation** Use CellGuide marker database → /opt/latch/plots-faas/runtime/mount/agent_config/context/technology_docs/xenium/cell_type_annotation.md
5. **Neighbors Enrichment Analysis** - Build spatial neighborhood graph → calculate centrality scores → calculate neighborhood enrichment score → identify spatially variable genes with **squidpy**.
6. **Domain Detection** — Convert the final `AnnData` object to `.h5ad` and launch `w_workflow(wf_name="wf.__init__.", ...)`.
7. **Cell Segmentation** — Ask the user whether to resegment. If confirmed, launch the workflow form by running `w_workflow(wf_name="wf.__init__.xenium_cell_segmentation_workflow", ...)`.

---

## Mandatory Package Redemption

Invoke the `redeem_package` tool to install required Xenium tools into the workspace.
  - `package_code`: `7a4f4bd980b3739a825072a975dd9a376c267ff7c84c1c9c59c8da196e58c3bd`
  - `package_version_id`: `401`

---

## General Workflow Rules

- **Always render a form** for any workflow with `lplots.widgets`.
  - Pre-populate sensible defaults.
  - **End with `w.value`** to execute.
- **Summarize results** (and progress/status) using `w_text_output(content=..., appearance={"message_box": “success” | “info”})`.
- **Render figures** using `w_plot` (and table widgets where appropriate). **DO NOT reuse the variable name fig. Use descriptive names: fig_qc, fig_umap, etc.**
  - Non-Scanpy figures: use **Plotly**.
- Do **not** ask users for file paths or node IDs; auto-discover from the attached folder.
- Keep objects and outputs reproducible: never overwrite without user consent; create new keys/versions.
- **Always include `w.value` at the end** when using a widget.
- Do **not** use `lplots.widgets.number`
- **ALWAYS** Wait for the workflow to complete before proceeding. After completion, always load the new adata from {output_directory}/{run_name} from ldata and never skip this.

---

## Step 0 — Data Transform

**Assumption (when this step is run):** The attached Xenium output directory contains the following files:

- `morphology_focus.ome.tif`
- `morphology_mip.ome.tif`
- `analysis.tar.gz`
- `transcripts.parquet`
- `cell_boundaries.parquet`
- `cell_feature_matrix.h5`
- `cells.parquet`

```python
params = {
    "input_file": LatchDir("latch://38438.account/Scratch/xenium/Input/Xenium_V1_FFPE_TgCRND8_17_9_months_outs"),
    "run_name": "my_run",
    "output_directory": LatchDir("latch:///Xenium_Preprocessing"),
}

w = w_workflow(
    wf_name="wf.__init__.xenium_preprocess_workflow",
    version=None,
    label="Launch Data Preparation Workflow",
    params=params,
    automatic=True,
)
execution = w.value
if execution is not None:
  res = await execution.wait()

  if res is not None and res.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
      # inspect workflow outputs for downstream analysis
      workflow_outputs = list(res.output.values())
```

---

## Step 1 — Data Loading

**Assumption:** Each sample already has a single, ready-to-use `.h5ad` containing counts, metadata, and spatial coordinates (either adata.obsm[Spatial] or adata.obsm[X_spatial]) generated from Step 0.

- Load the .h5ad file and render it with w_h5(adata, spatial_dir), where spatial_dir is the directory containing the .pmtiles files (should be the same directory as the .h5ad).
- If embeddings or clusters already exist in adata, do not run preprocessing yet. Ask for confirmation in **Step 2** before recomputing.
- If multiple samples are provided
  - Merge the .h5ad files into a single AnnData object and add adata.obs["batch_key"] to preserve sample information.
  - Create adata.obsm["Offset"] as a shifted spatial embedding where each sample’s coordinates are translated by a sample-specific offset so they appear side by side in a shared space instead of overlapping.

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
- **Batch integration**
  - If multiple samples are provided, integrate them using Scanpy `sce.pp.harmony_integrate(key="batch_key")`

- **Neighborhood graph, UMAP, clustering**
  - Build neighbors on the integrated space.
  - Compute UMAP; store in `adata.obsm["X_umap"]` or `adata.obsm["UMAP"]`.
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
- **Default custering column:** `K=5`. Offer a dropdown widget to let user select other number of clusters.
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
- Select the top four biologically meaningful marker genes and color the spatial embedding by their log1p expression in four subplots. Explain why they are biologically meaningful.

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

- Use DGE markers and the CellGuide per-gene database to assign a **biological** cell type to each cluster.
- If cluster values are numerical, convert to string before assigning cell types.
- Use CellGuide per-gene database to propose raw cell types per cluster, then map them into a **controlled vocabulary**:
  - If `/opt/latch/plots-faas/runtime/mount/agent_config/context/technology_docs/xenium/cell_type_vocab_index.json` is present and a matching vocab config exists for the dataset’s organism and tissue, treat that config’s `allowed_vocab` as the **only** valid labels and map CellGuide types to these labels.
  - If no matching vocab config is found, use cleaned CellGuide cell-type names directly as labels.
- **NEVER** use placeholder or vague labels such as `cluster_0`, `Cell_Population_1`, `celltype_placeholder_0`, `label_1`, `unknown`, `other`, `neuron associated cell`.
- For each cluster, store the final label, confidence, key markers and supporting CellGuide types in `adata.uns["cell_type_annotation"]`, and write labels to `adata.obs["cell_type"]` as a categorical, updating the existing `w_h5` viewer in place.
- Visualize annotated cell types in both the spatial embedding and the UMAP embedding subplots.
- Verify that `adata.obs["cell_type"]` contains no missing values and the labels make sense. If it is missing or contains nulls, correct this before proceeding.

---

## Step 5 — Neighbors Enrichment Analysis

- Build spatial neighborhoold graph using clustering results and calculate centrality score. Give a one-sentence summary of what each score represents and how to interpret it with `w_text_output`
- Visualize centrality scores in spatial embeddings.
```python
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True, spatial_key="Spatial")
sq.gr.centrality_scores(adata, cluster_key="leiden")
sq.pl.centrality_scores(adata, cluster_key="leiden", figsize=(20, 5))

fig_centrality = plt.gcf()
w_plot(label="Centrality Scores by Cell Type", source=fig_centrality)
plt.close()
```
- Calculate enrichment score using `sq.pl.nhood_enrichment` and visualize the z-score in Plotly heatmap.
```python
enrichment_zscore = adata.uns['cell_type_nhood_enrichment']['zscore']
cell_type_names = adata.obs['cell_type'].cat.categories.tolist()

# Create heatmap with plotly
fig_enrichment = go.Figure(data=go.Heatmap(
    z=enrichment_zscore,
    x=cell_type_names,
    y=cell_type_names,
    colorscale='RdBu_r',
    zmid=0,
    colorbar=dict(title="Z-score"),
    text=enrichment_zscore.round(2),
    texttemplate='%{text}',
    textfont={"size": 10}
))

fig_enrichment.update_layout(
    title="Neighborhood Enrichment Analysis",
    xaxis_title="Cell Type",
    yaxis_title="Cell Type",
    width=700,
    height=650
)

w_plot(label="Neighborhood Enrichment Heatmap", source=fig_enrichment)
```
- Identify spatially variable genes by Moran's I score and and color the spatial embedding by top four genes' expression in four subplots.
```python
sq.gr.spatial_autocorr(
    adata,
    mode="moran",
    n_jobs=-1,
)
```
---

## Step 6 — Domain Detection (Optional)

- Rename adata.obsm["Spatial"] to adata.obsm["X_spatial"], then save the processed object to `.h5ad`.
    ```python
    from latch.ldata.path import LPath
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
- Load {output_dir}/*.h5ad from ldata to `adata_domain_detected` and colour cells by `adata_domain_detected.obs[labels_scaled_gaussian_*]` in spatial embedding.

```python
from latch.types import LatchDir
params = {
    "input_file": LatchFile("latch:///xenium_processed_for_domain_detection.h5ad"),
    "run_name": "my_run",
    "output_dir": LatchDir("latch:///Domain_detection_output"),
    "lambda_list": "0.5",
}

w = w_workflow(
    wf_name="wf.__init__.domain_detection_wf",
    version=None,
    label="Launch Domain Detection Workflow",
    params=params,
    automatic=True,
)
execution = w.value

if execution is not None:
  res = await execution.wait()

  if res is not None and res.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
      # inspect workflow outputs for downstream analysis
      workflow_outputs = list(res.output.values())
```
---

## Step 7 — Cell Segmentation (Optional)

- Ask: “Re-segment cells?”
- If **yes**:
  - **Use the full-resolution morphology image only**: `morphology.ome.tif` (never any file containing “mip”).
  - **Render a form** and launch
    `w_workflow(wf_name="wf.__init__.xenium_cell_segmentation_workflow", …).value`
  - **Expose parameters as widgets** (at minimum: `output_dir`, `level`) with sensible defaults.
- If **no**, proceed without changes.
- Wait for the workflow to complete.

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
    label="Launch Cell Segmentation Workflow",
    params=params,
    automatic=True,
)
execution = w.value
```

---
