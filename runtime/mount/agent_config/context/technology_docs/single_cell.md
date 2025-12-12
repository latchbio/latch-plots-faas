# 10x Genomics scRNA-seq Analysis Guide

---

<analysis_overview>

## Pipeline Overview

0. **Data Loading**
1. **Quality Control (QC)**
2. **Preprocessing**
3. **Clustering and Differential Gene Expression Analysis**
4. **Cell Type Annotation**
5. **Gene Set Enrichment Analysis**
6. **Summary and Explorer**

</analysis_overview>



---



<data_loading>

## Step 0 - Data Loading

**Goals:** Create an analysis header, and read the 10x Genomics h5 data format as an AnnData using ScanPy.

```python
sc.read_10x_h5()
```

</data_loading>



<quality_control>

## Step 1 - Quality Control (QC)

**Goal:** Create a new tab and analysis header. Run QC metrics on the AnnData and display the results through visualizations.


First identify mitochondrial genes, ribosomal genes, and hemoglobin genes.
```python
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
```

Then, compute standard QC metrics with ScanPy.
```python
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
```

Start by displaying the ```n_genes_by_counts```, ```total_counts```, and ```pct_counts_mt``` to the user as violin plots through a widget with ```w_plot```.

```python
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)
```

**Render QC parameters as widgets with sensible defaults**: ```min_genes``` (default: 100), ```min_cells``` (default: 3)

> Generate widgets for the parameters above; pre-fill with defaults. Require explicit user confirmation before running.

Before moving on, explicitly ask the user if they would like to further explore the quality of the single-cell data, including:
- Visualizing percent ribosomal counts or percent hemoglobin counts through violin plots (```pct_counts_ribo``` or ```pct_counts_hb```)
- Scatter plots of joint QC variables (ex. ```sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")```)

</quality_control>



<preprocess_data>

## Step 2 - Preprocessing

**Goal:** Create a new tab and analysis header. Perform normalization and dimensionality reduction

- Normalization
    - Do library-size normalization, then log transform.
- Dimensionality reduction
    - PCA on all genes, store in ```adata.obsm["X_pca"]```
- Neighborhood graph, UMAP, clustering
    - Build neighbors on the integrated space.
    - Compute UMAP; store in adata.obsm["X_umap"] or adata.obsm["UMAP"].

**Render parameters as widgets with sensible defaults.**: 
- **PCA**: ```n_pcs``` (default: 50)
- **Neighbors:** ```n_neighbors``` (default: 15), ```metric``` (default: ```"euclidean"```)

> Generate widgets for the parameters above; pre-fill with defaults. Require explicit user confirmation before running.

</preprocess_data>



<clustering_and_dge>

## Step 3 - Clustering and Differential Gene Expression

**Goal:** Create a new tab and analysis header. Run Leiden clustering and perform differential gene expression on each cluster.

- **Leiden Clustering**
    - Run Leiden clustering; store categorical labels in adata.obs["leiden"].
    - Render Leiden clusters on the spatial embedding and UMAP using Plotly subplots.
- **Differential Gene Expression (DGE)**
    - Run DGE with ScanPy:
```sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')```
    - Report top marker genes for each cluster and create a dot plots with ScanPy.
    - Select the top four biologically meaningful marker genes and color the spatial embedding by their log1p expression. Explain why they are biologically meaningful.

</clustering_and_dge>



<cell_type_annotation>

## Step 4 - Cell Type Annotation

**Goal:** Create a new tab and analysis header. Annotate the cell types of the AnnData.

- Use the DGE markers and the CellGuide per-gene database to assign a **biological** cell type to each cluster.
- If cluster values are numerical, convert to string before assigning cell types.
- Use CellGuide per-gene database to propose raw cell types per cluster.
- **NEVER** use placeholder or vague labels such as `cluster_0`, `Cell_Population_1`, `celltype_placeholder_0`, `label_1`, `unknown`, `other`, `neuron associated cell`.
- For each cluster, store the final label, confidence, key markers and supporting CellGuide types in `adata.uns["cell_type_annotation"]`, and write labels to `adata.obs["cell_type"]`.
- Visualize annotated cell types in the UMAP embedding plot.
- Verify that `adata.obs["cell_type"]` contains no missing values and the labels make sense. If it is missing or contains nulls, correct this before proceeding.

</cell_type_annotation>



<gsea_analysis>
## Step 5 - Gene Set Enrichment Analysis

**Goal:** Create a new tab and analysis header. Identify pathways that are highly enriched in each ```leiden``` cluster and explain the biological significance.

- Use GSEApy to do gene set enrichment for each cluster. Use GSEApy gene sets that are biologically relevant to the given dataset.

> Use Enrichr API directly instead of the gseapy module, which isn’t available

- **Create a dotplot of the top enriched terms per cluster.**



</gsea_analysis>



<summary_and_exploration_analysis>

## Step 6 - Summary and Exploration

**Goal:** Create a new tab and analysis header. Provide an overview of analysis steps and display an interactive AnnData viewer with ```w_h5```

- Create a table with a complete pipeline overview along with key outputs.
- Create a second table to summarize the cell type annotations with four columns:
    - Cell Type
    - Count
    - Percent
    - Key Markers
- **Create interactive AnnData explorer using ```w_h5``` viewer and UMAP embeddings**. The user can:
    - Color by cell type or leiden cluster
    - View gene expression for specific genes
    - Lasso select cells for further analysis

</summary_and_exploration_analysis>
