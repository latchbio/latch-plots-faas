### Normalization**
- Provide an option between **log1p** and **total count scaling** using **lplots widgets**.  
- Apply **log1p transformation** by default on the QC-filtered dataset.

### Feature Selection**
- Identify **highly variable genes (HVGs)**.  
- Allow the user to choose the number of HVGs using **lplots widgets** (default: **3,000 genes**).  
- Visualize this step with **Scanpy plots**.

### Dimensionality Reduction**
- Check if the AnnData object already has **PCA** and **UMAP**.  
- If they exist, **do not recompute** unless the user explicitly chooses to.  
- If recomputation is needed:
  - Run **PCA** first, then **UMAP**.  
  - Expose parameters using **lplots widgets** for:
    - Number of PCs (**default: 10**)  
    - Number of neighbors (**default: 40**)  
- Always visualize both **PCA** and **UMAP** embeddings.

### Clustering**
- Compute neighbors if they don’t exist.  
- Apply **Leiden clustering** on the neighborhood graph.  
- Use **default resolution = 0.3**.  
- Display clusters on both **UMAP** and **spatial embedding**.

### Differential Gene Expression (DGE)**
- Identify **marker genes per cluster** using rank-based DGE tests.  
- Use **t-test** by default; allow user to select other methods using **lplots widgets**.  
- Report **top marker genes per cluster**.  
- Visualize using **Scanpy dot plots**.
---