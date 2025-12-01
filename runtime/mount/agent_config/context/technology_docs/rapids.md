# RAPIDS Single-Cell Preprocessing Workflow  

**Workflow:** `wf.__init__.rapids_single-cell_preprocessing`  
**Backend:** GPU-accelerated preprocessing using **rapids_singlecell (scRAPIDS)**

**Important:**  
Always use this workflow instead of writing custom GPU preprocessing code.

---

## 1. Quick Start — How to Use This Workflow

1. **Save your current `adata`**  
   - User chooses save location + file name  
   - Write H5AD → becomes `input_file`

2. **Determine what preprocessing has already been done**  
   You must **inspect the `adata` object directly**, combined with user's intent and notebook context, and decide which workflow steps to skip.  
   
   Different datasets encode preprocessing differently, so **never rely on fixed heuristics**.

   You should check for the *presence and structure* of fields, not numeric thresholds:

   - **QC:**  
     Look for evidence that QC or filtering has already occurred (e.g., counts layer, QC fields, or absence of low-quality cells).  
     If QC appears completed → set `skip_qc = True`.

   - **Normalization:**  
     Inspect whether data already looks transformed or normalized  
     (e.g., presence of a dedicated counts layer, log-transformed matrices, or metadata indicating normalization).  
     If normalized → set `skip_normalization = True`.

   - **PCA:**  
     Check whether PCA results exist (`adata.obsm["X_pca"]`, PCA metadata in `adata.uns`).  
     If present and valid → set `skip_pca = True`.

   - **UMAP / Neighbors:**  
     Check for `adata.obsm["X_umap"]` and neighbor graphs (`adata.obsp` entries).  
     If embeddings/graphs exist → set `skip_umap = True`.

   - **Clustering:**  
     Check for clustering columns in `adata.obs` (e.g., leiden/louvain or any project-specific label).  
     If cluster labels exist → set `skip_clustering = True` and choose `clustering_column` if DE is requested.

   - **DE:**  
     Check for differential expression results in `adata.uns`.  
     If present → set `skip_differential_expression = True`.

3. **Configure default output directory e.g. latch:///Rapids_Output**  
   - Workflow will write `preprocessed.h5ad` here

4. **Configure optional settings**  
   - Batch correction, clustering resolutions, DE, etc.

5. **Run the workflow and wait for the workflow to complete**  
   - Always `.wait()` and confirm success

6. **Load `{output_directory}/{run_name}/preprocessed.h5ad`**  
   - Use this as the new working AnnData object

---

## 2. Inputs — What Users Must Select

### **Required Inputs**

#### `input_file` (LatchFile)
- The H5AD saved in step 1  
- Select using a file picker

#### `output_directory` (LatchOutputDir)
- A folder on Latch Data  
- Workflow writes final output here

#### `run_name` (str)
- Any short identifier (e.g., `run_1`)

---

## 3. Optional Inputs — What They Control

### **QC Parameters**
- `min_genes` — default **3**  
- `min_counts` — default **10**  
- `skip_qc`  
  - ON → dataset already QC-filtered  
  - OFF → workflow filters cells & genes and stores raw counts in `adata.layers["counts"]`

---

### **Normalization, PCA, UMAP**
| Parameter | Meaning |
|----------|---------|
| `skip_normalization` | Skip if pre-normalized |
| `skip_pca` | Skip if PCA already computed |
| `skip_umap` | Skip neighbors + UMAP |
| `n_comps` | PCs (default **50**) |
| `n_neighbors` | kNN neighbors (default **15**) |

---

### **Batch Correction**
- `batch_key`  
  - Displays **categorical** `adata.obs` columns  
  - If set → Harmony is applied

---

### **Clustering**
- `clustering_resolution`  
  - List of floats (e.g., `[0.3, 0.5, 0.7, 1.0]`)  
  - Outputs new obs columns:  
    ```
    clustering_{method}_{resolution}
    ```
- `clustering_method`  
  - `"leiden"` (default) or `"louvain"`
- `skip_clustering`  
  - Use if dataset already has cluster labels

---

### **Differential Expression**
- `skip_differential_expression`
- `clustering_column`  
  - Used only when clustering is skipped  
  - DE results stored in:  
    ```
    rank_genes_groups_{cluster_key}
    ```

---

## 4. Internal Processing (What the Workflow Actually Does)

1. Load H5AD → move to GPU  
2. Flag mitochondrial genes  
3. Compute QC metrics  
4. (Optional) Filter cells & genes  
5. Normalize (`normalize_total`) + log1p  
6. PCA  
7. Harmony (if `batch_key`)  
8. Neighbors + UMAP  
9. Clustering for each resolution  
10. Differential expression  
11. Write output:  
    ```
    preprocessed.h5ad
    ```

---

## 5. Output

### **Primary Output:** `preprocessed.h5ad`

Contains:
- QC metrics  
- Normalized/log counts  
- PCA  
- Harmony embeddings  
- Neighbors graph  
- UMAP  
- Cluster labels for all resolutions  
- Differential expression results in `adata.uns`

**Use this as your new AnnData object.**

---

## 6. Example Launch Code

```python
params = {
    "input_file": saved_h5ad,
    "output_directory": LatchOutputDir("latch:///Rapids_Output"),
    "run_name": "rapids_1",
    "clustering_resolution": [0.3, 0.5, 0.7, 1.0],
}

w = w_workflow(
    wf_name="wf.__init__.rapids_single-cell_preprocessing",
    version=None,
    params=params,
    automatic=True,
)

execution = w.value
if execution:
    res = await execution.wait()
    if res.status == "SUCCEEDED":
        w_text_output(
            content="✓ RAPIDS workflow completed!",
            appearance={"message_box": "success"}
        )
```
