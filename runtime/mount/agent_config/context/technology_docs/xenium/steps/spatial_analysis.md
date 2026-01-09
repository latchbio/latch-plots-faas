# Step 5 — Neighbors Enrichment Analysis and Domain Detection

<goal>
Quantify spatial neighborhood structure (centrality, enrichment, spatially variable genes) and detect higher-order spatial domains from the processed Xenium dataset.
</goal>

<method>
1/ **Neighbors Enrichment and Spatial Structure**: Build a spatial neighborhood graph using:
   - `sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True, spatial_key="Spatial")`
   Then compute **centrality scores** with `sq.gr.centrality_scores` and visualize them on the spatial embedding, providing a short explanation of what each centrality score represents and how to interpret it.

2/ Compute neighborhood enrichment scores between cell types with:
   - `sq.gr.nhood_enrichment(adata, cluster_key="cell_type")`
   Visualize the resulting z-scores in a Plotly heatmap.

3/ Identify spatially variable genes using Moran’s I with:
   - `sq.gr.spatial_autocorr(adata, mode="moran", n_jobs=-1)`
   Highlight the top spatially autocorrelated genes on the spatial embedding.

4/ **Domain Detection**: Run the domain detection workflow as specified in `wf/domain_detection_wf.md` using the preprocessed `.h5ad` as input, and load the resulting domain-labeled `.h5ad` for downstream visualization.
</method>

<workflows>
- `wf/domain_detection_wf.md`
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
- The processed `.h5ad` is successfully written, uploaded, and consumed by the domain-detection workflow.
- Domain labels (e.g. `labels_scaled_gaussian_*`) are present in the output `AnnData` and produce coherent spatial domains when visualized.
</self_eval_criteria>
