# Step 3 — Differential Gene Expression (DGE)

<goal>
Identify cluster-specific marker genes suitable for annotation and visualization.

</goal>

<method>

- **Default test (fastest):** `t-test_overestim_var` with Benjamini–Hochberg FDR control.
- **Alternative (more robust to outliers/zero inflation):** `wilcoxon`
- **Default clustering column:** `"K=5"`. Offer a dropdown widget to let the user select other cluster columns.

> **Heuristic:** keep `t-test_overestim_var` as default (“Quick mode”). Offer a toggle to “Robust mode” → `wilcoxon` when users prioritize robustness over speed.

### Parameters (render as widgets with defaults)
- **Method:** `t-test_overestim_var` *(default)* \| `wilcoxon`
- **Top genes per cluster:** `top_n` *(default: 5)*

### Core DGE call

```python
valid_cells = ~adata.obs[cluster_col].isna()
n_removed = (~valid_cells).sum()

adata_filtered = adata[valid_cells].copy()
adata_filtered.obs[cluster_col] = adata_filtered.obs[cluster_col].apply(
    lambda x: str(x) if not isinstance(x, str) else x
)

sc.tl.rank_genes_groups(
    adata_filtered,
    groupby=cluster_col,
    method=dge_method,
    use_raw=False,
    corr_method="benjamini-hochberg",
)
```

### Reporting

- Per cluster, report:
    - gene
    - log fold change
    - p-value
    - **FDR** (adjusted p-value).
- Store results under `adata.uns["dge"]`.
- Report top marker genes for each cluster and make dot plots with Scanpy.
- Select the top four biologically meaningful marker genes and color the spatial embedding by their `log1p` expression in four subplots. Explain briefly why they are biologically meaningful.

```python
sc.pl.dotplot(
    adata,
    var_names=all_marker_genes,
    groupby=cluster_col_name,
    standard_scale="var",
    show=False,
)

dotplot_fig = plt.gcf()
w_plot(label="Marker Gene Expression Dot Plot", source=dotplot_fig)
plt.close()
```
</method>

<workflows>
</workflows>

<library>
- `scanpy`
- `matplotlib`
- `plotly`
</library>

<self_eval_criteria>

-  DGE runs without errors on the selected cluster column.
- Most clusters have non-empty marker sets with reasonable log fold changes and FDR.
- `adata.uns["dge"]` is populated and can be reused for cell-type annotation.
- Spatial and UMAP plots of top markers look biologically plausible.
</self_eval_criteria>

