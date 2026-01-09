# Step 5 — Neighbors Enrichment Analysis and Domain Detection

<goal>
Quantify spatial neighborhood structure (centrality, enrichment, spatially variable genes) and detect higher-order spatial domains from the processed Xenium dataset.
</goal>

<method>
## Neighbors Enrichment and Spatial Structure

- Build a spatial neighborhood graph using clustering results and compute centrality scores.
- Use `w_text_output` to provide a one-sentence summary of what each centrality score represents and how to interpret it.
- Visualize centrality scores on spatial embeddings.

```python
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True, spatial_key="Spatial")
sq.gr.centrality_scores(adata, cluster_key="leiden")
sq.pl.centrality_scores(adata, cluster_key="leiden", figsize=(20, 5))

fig_centrality = plt.gcf()
w_plot(label="Centrality Scores by Cell Type", source=fig_centrality)
plt.close()
````

- Calculate neighborhood enrichment scores with Squidpy and visualize the z-scores in a Plotly heatmap.

```python
sq.gr.nhood_enrichment(adata, cluster_key="cell_type")
enrichment_zscore = adata.uns["cell_type_nhood_enrichment"]["zscore"]
cell_type_names = adata.obs["cell_type"].cat.categories.tolist()

fig_enrichment = go.Figure(
    data=go.Heatmap(
        z=enrichment_zscore,
        x=cell_type_names,
        y=cell_type_names,
        colorscale="RdBu_r",
        zmid=0,
        colorbar=dict(title="Z-score"),
        text=enrichment_zscore.round(2),
        texttemplate="%{text}",
        textfont={"size": 10},
    )
)

fig_enrichment.update_layout(
    title="Neighborhood Enrichment Analysis",
    xaxis_title="Cell Type",
    yaxis_title="Cell Type",
    width=700,
    height=650,
)

w_plot(label="Neighborhood Enrichment Heatmap", source=fig_enrichment)
```

* Identify spatially variable genes using Moran's I and color the spatial embedding by the top four genes’ expression in four subplots.

```python
sq.gr.spatial_autocorr(
    adata,
    mode="moran",
    n_jobs=-1,
)
```

## Domain Detection

- Ensure spatial coordinates are stored in `adata.obsm["X_spatial"]` (renaming from `"Spatial"` if needed).
- Save the processed `AnnData` object to an `.h5ad` file and upload it to Latch.

```python
from pathlib import Path
from latch.ldata.path import LPath

if "Spatial" in adata.obsm:
    adata.obsm["X_spatial"] = adata.obsm["Spatial"]

output_h5ad_path = "/tmp/xenium_processed_for_domain_detection.h5ad"
adata.write_h5ad(output_h5ad_path)

output_lpath = LPath("latch:///xenium/analysis_output/xenium_processed_for_domain_detection.h5ad")
local_path = Path(output_h5ad_path)
output_lpath.upload_from(local_path)
```

- Use this `.h5ad` as input to the domain-detection workflow and wait for completion.

```python
from latch.types import LatchDir, LatchFile
from lplots.widgets.workflow import w_workflow

params = {
    "input_file": LatchFile("latch:///xenium/analysis_output/xenium_processed_for_domain_detection.h5ad"),
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
        workflow_outputs = list(res.output.values())
```

- Load `{output_dir}/*.h5ad` from `ldata` into `adata_domain_detected` and color cells in the spatial embedding by `adata_domain_detected.obs["labels_scaled_gaussian_*"]` (or equivalent domain label columns).

</method>

<workflows>
- `wf.__init__.domain_detection_wf`
</workflows>

<library>
- `squidpy`
- `plotly`
- `matplotlib`
- `scanpy`
- `latch.types`
- `latch.ldata.path`
- `lplots.widgets.workflow`
- `pathlib`
</library>

<self_eval_criteria>

- Spatial neighbor graph builds successfully with reasonable centrality patterns.
- Neighborhood enrichment heatmap shows interpretable positive and negative z-scores between cell-type pairs.
- Moran’s I identifies non-trivial sets of spatially variable genes.
- The processed `.h5ad` is successfully written, uploaded, and accepted by the domain-detection workflow.
- Domain labels (e.g. `labels_scaled_gaussian_*`) are present in the output `AnnData` and produce coherent spatial domains when visualized.
</self_eval_criteria>

```
