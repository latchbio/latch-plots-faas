<method>

**Quality control**
- Compute standard QC metrics (genes per cell, counts per cell, mitochondrial fraction if available).
```python
sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
```
- Show histograms using plotly; expose filter thresholds via widgets `w_text_input`.
- Confirm before applying filters; keep an unfiltered copy.
- **ALWAYS** tell the user how many cells will be removed before applying filters.
- If all cells are removed, ask user to re-enter parameters.

**Normalization**
- Library-size normalization, then log transform.

**Dimensionality reduction**
- PCA on all genes; store in `adata.obsm["X_pca"]`.
```python
sc.tl.pca(adata, n_comps=n_pcs)
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=n_pcs, show=False)

# Capture the figure and render with the plot widget
fig_pca = plt.gcf()
w_plot(label="PCA Explained Variance Ratio", source=fig_pca)

plt.close()
```
**Batch integration**
- If multiple samples are provided, integrate them using Scanpy `sce.pp.harmony_integrate(key="batch_key")`

**Neighborhood graph, UMAP, clustering**
- Build neighbors on the integrated space.
- Compute UMAP; store in `adata.obsm["X_umap"]` or `adata.obsm["UMAP"]`.
- Run Leiden clustering; store categorical labels in adata.obs["leiden"]; **update the original w_h5 viewer in place—do not create a new viewer.**
- Render Leiden clusters on the spatial embedding and UMAP using Plotly subplots.

**Non-destructive policy**
- Do not overwrite existing embeddings or clusters unless the user opts in. Write to new keys if needed.

### Parameters (render as widgets with sensible defaults)
- **QC thresholds:** `min_genes` (default: 5), `min_counts` (default: 20), `max_mt_frac` (default: 0.2)
- **PCA:** `n_pcs` (default: 50)
- **Neighbors:** `n_neighbors` (default: 15), `metric` (default: `"euclidean"`)
- **Leiden:** `resolution` (default: 0.5), `random_state` (default: 0)

> Generate widgets for the parameters above; pre-fill with defaults. Require explicit user confirmation before running.

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
</method>
