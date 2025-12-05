# Cell2location Spatial Transcriptomics Deconvolution Workflow  

**Workflow:** `wf.__init__.cell2location`  
**Backend:** Bayesian deconvolution using **cell2location** to map cell types from single-cell reference onto spatial transcriptomics data

**Important:**  
Always use this workflow instead of writing custom cell2location code for spatial transcriptomics deconvolution.

---

## 1. Quick Start — How to Use This Workflow

1. **Prepare your spatial transcriptomics data (`query_adata`)**  
   - Save spatial transcriptomics AnnData to Latch Data → becomes input to `query_adata`
   - **Requirement:** Must contain **raw counts in `adata.X`** for compatibility

2. **Prepare your single-cell reference data (`ref_adata`)**  
   - **If user has an existing H5AD file:**
     - Ask user to provide the H5AD file
     - **Verify** that it contains:
       - **Raw counts in `adata.X`** 
       - **Cell type labels in `adata.obs`** 
     - If verification passes → use this file as `ref_adata`
     - If verification fails → inform user, ask for alternative and offer to help curate from a paper. 
   - **If user does not have a reference:**
     - Curate a single-cell reference from scratch following the workflow in `technology_docs/tools/curate_sc_reference.md`
     - Save the curated reference → becomes `ref_adata`

3. **Determine batch keys (if applicable)**  
   - Check if your datasets have batch effects that need correction
   - Inspect `query_adata.obs` and `ref_adata.obs` for batch columns
   - If batches exist → set `query_batch_key` and/or `ref_batch_key`

4. **Identify the cell type column in reference**  
   - Inspect `ref_adata.obs` to find the column containing cell type labels
   - Common names: `"cell_type"`, `"celltype"`, `"annotation"`, etc.
   - Set `labels_key` to this column name

5. **Configure output directory**  
   - Default: `latch:///Cell2location`
   - Workflow will write deconvolution results here

6. **Configure training and inference parameters**  
   - Adjust `max_epochs`, `batch_size`, `num_samples` based on dataset size
   - Set `N_cells_per_location` based on expected cell density
   - Configure `detection_alpha` for detection sensitivity

7. **Run the workflow and wait for completion**  
   - Always `.wait()` and confirm success

8. **Load `{output_directory}/{run_name}/sp_query.h5ad`**  
   - Contains spatial data with cell type proportions per spot
   - Use this as the new working AnnData object

---

## 2. Inputs

### **Required Inputs**

#### `query_adata` (LatchFile)
- The spatial transcriptomics H5AD 
- **Must contain raw counts in `adata.X`**

#### `ref_adata` (LatchFile)
- The single-cell reference H5AD
- **Must contain raw counts in `adata.X` and cell type labels in `adata.obs`**

#### `output_directory` (LatchOutputDir)
- A folder on Latch Data  
- Workflow writes final output here
- Default: `latch:///Cell2location`

#### `run_name` (str)
- Any short identifier (e.g., `bone_visium`)

#### `labels_key` (str)
- Name of the `ref_adata.obs` column containing cell type labels

---

## 3. Optional Inputs — What They Control

### **Batch Correction**
- `query_batch_key`  
  - Name of categorical column in `query_adata.obs` that references samples, conditions, or batch. 
  - Set if your spatial data contains multiple samples, batches, or experimental conditions
  - Leave as `None` if spatial data comes from a single sample/condition

- `ref_batch_key`  
  - Name of categorical column in `ref_adata.obs` that references samples, conditions, or batch. 
  - Set if your reference data contains multiple samples, batches, or experimental conditions
  - If set → batch effects in reference data are corrected during model training
  - Leave as `None` if reference data comes from a single sample/condition

---

### **Training Parameters**
| Parameter | Meaning | Default |
|----------|---------|---------|
| `max_epochs` | Maximum training epochs | **300** |
| `batch_size` | Batch size for training | **2500** |
| `num_samples` | Number of posterior samples | **1000** |
| `train_size` | Fraction of data for training (0-1) | **1.0** (use all data) |
| `n_top_genes` | Number of highly variable genes to use | **2000** |

---

### **Inference Parameters**
- `N_cells_per_location`  
  - Expected number of cells per spatial location/spot  
  - Default: **15**  
  - Adjust based on your technology:
    - Visium (10x): ~10-20 cells per spot
    - Slide-seq: ~5-15 cells per spot
    - MERFISH: varies by resolution

- `detection_alpha`  
  - Detection sensitivity threshold  
  - Default: **20**  
  - Lower values = more sensitive (detect more cell types)
  - Higher values = more conservative (only strong signals)

---

## 4. Outputs

The workflow writes output to `latch:///{output_directory}/{run_name}/` with the following structure:

### **Primary Output Files:**

#### `sp_query.h5ad`
- Spatial transcriptomics data with deconvolution results

#### `sc_reference.h5ad`
- Processed reference single-cell data used for training
- May contain additional annotations added during model training

#### `model.pt`
- Trained cell2location model weights
- Can be used for inference on additional spatial datasets without retraining

### **Output Directories:**

#### `query_model/`
- Model artifacts and intermediate files for the query (spatial) data

#### `reference_model/`
- Model artifacts and intermediate files for the reference (single-cell) data

### **Logs:**

#### `task_log.txt`
- Workflow execution log
- Contains training progress, warnings, and error messages

---

## 6. Example Launch Code

```python
spatial_h5ad = LatchFile("latch://<account_id>/spatial.h5ad")
reference_h5ad = LatchFile("latch://<account_id>/sc_ref.h5ad")

params = {
    "query_adata": spatial_h5ad,
    "ref_adata": reference_h5ad,
    "output_directory": LatchOutputDir("latch:///Cell2location"),
    "run_name": "bone_visium",
    "labels_key": "cell_type",
    "N_cells_per_location": 15,
    "n_top_genes": 2000,
    "max_epochs": 300,
    "batch_size": 2500,
    "num_samples": 1000,
    "detection_alpha": 20,
    "query_batch_key": None,  # or "batch" if batch correction needed
    "ref_batch_key": None,    # or "batch" if batch correction needed
    "train_size": 1.0,
}

w = w_workflow(
    wf_name="wf.__init__.cell2location",
    version=None,
    params=params,
    automatic=True,
)

execution = w.value
if execution:
    res = await execution.wait()
    if res.status == "SUCCEEDED":
        w_text_output(
            content="✓ Cell2location workflow completed!",
            appearance={"message_box": "success"}
        )
```
