<goal>
Perform dimensionality reduction using PCA and UMAP, and apply Leiden clustering on normalized Pixelator PNA data.
</goal>

<method>
1/ **Parameter Selection**
   - Allow user to specify:
     - Number of top highly variable genes (default: "50")
     - Number of PCs for UMAP calculation (default: "10")
     - Number of neighbors for UMAP calculation (default: "30")
     - Leiden clustering resolution (default: user-specified)

2/ **Highly Variable Genes and Scaling**
   - Import required libraries:
     ```python
     import scanpy as sc
     import numpy as np
     ```
   - Calculate highly variable genes:
     ```python
     sc.pp.highly_variable_genes(
         adata, flavor="seurat_v3", n_top_genes=int(top_genes)
     )
     ```
   - Scale the DSB layer:
     ```python
     adata.layers["scaled_dsb"] = sc.pp.scale(
         adata, zero_center=True, layer="dsb", copy=True
     ).layers["dsb"]
     ```

3/ **PCA Computation**
   - Compute PCA:
     ```python
     adata.obsm["pca"] = sc.tl.pca(
         adata.layers["scaled_dsb"], random_state=42
     )
     ```

4/ **Neighbors Graph Construction**
   - Compute neighbors:
     ```python
     sc.pp.neighbors(
         adata,
         n_neighbors=int(neighbors),
         n_pcs=int(npcs),
         use_rep="pca",
         metric="cosine",
     )
     ```

5/ **UMAP Embedding**
   - Compute UMAP:
     ```python
     sc.tl.umap(
         adata,
         min_dist=0.3,
     )
     ```

6/ **Leiden Clustering**
   - Apply Leiden clustering:
     ```python
     sc.tl.leiden(adata, resolution=float(leiden_res))
     n_clusters = adata.obs['leiden'].nunique()
     ```

7/ **Visualization**
   - Display the adata object using an h5 viewer
</method>

<workflows>
</workflows>

<library>
</library>

<self_eval_criteria>
</self_eval_criteria>
