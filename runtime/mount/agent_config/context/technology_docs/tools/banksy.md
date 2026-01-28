# BANKSY Spatial Domain Detection Workflow

**Workflow:** `wf.__init__.domain_detection_wf`  
**Backend:** Graph-based spatial domain detection using **BANKSY**

**Purpose:**  
Detect spatial domains (tissue regions / niches) from **preprocessed spatial transcriptomics data**.  
This workflow is **technology-agnostic** and works with any spatial transcriptomics dataset that provides gene expression and spatial coordinates.

Always use this workflow instead of writing custom BANKSY code.

---

## 1. When to Use This Workflow

Use this workflow when you want to:
- Identify **spatial domains** from spatial transcriptomics data
- Assign **domain labels per observation** and store them in `adata.obs`
- Analyze tissue architecture, spatial niches, or regional programs

**Compatible technologies include (but are not limited to):**
- Xenium
- Visium
- Slide-seq / Slide-seqV2
- MERFISH
- CosMx
- seqFISH
- Any spatial dataset stored as a valid `.h5ad` with coordinates

**Prerequisite:**  
Input `.h5ad` must already be **preprocessed** (QC, normalization, filtering).

---

## 2. Required Inputs

### `input_file` (`LatchFile`)
- Preprocessed spatial transcriptomics `.h5ad`
- Must contain spatial coordinates in:
  - `adata.obsm["X_spatial"]`

If your data uses a different key (for example `Spatial`), map it before saving:
```python
if "X_spatial" not in adata.obsm and "Spatial" in adata.obsm:
    adata.obsm["X_spatial"] = adata.obsm["Spatial"]
```

---

### `run_name` (`str`)
- Short identifier for this run  
- Used to create `{output_dir}/{run_name}/`

---

### `output_dir` (`LatchDir`)
- Base output directory on Latch
- Domain detection results are written under:
```
{output_dir}/{run_name}/
```

---

### `lambda_list` (`str`)
- **Single lambda value encoded as a string**
- Controls spatial smoothness in BANKSY

Examples:
- `"0.3"`
- `"0.5"`
- `"0.8"`

> Note: Despite the parameter name, the workflow currently accepts **one lambda value per run**.

---

## 3. Outputs

Written under:
```
{output_dir}/{run_name}/
```

### Primary Outputs
- One `.h5ad` file per run
- Contains **domain labels in `adata.obs`**
  - Example column name:
    - `labels_scaled_gaussian_0.5`

The output AnnData is ready for downstream visualization and analysis.

---

## 4. Minimal Example — Launch Workflow

```python
from pathlib import Path
from latch.ldata.path import LPath
from latch.types import LatchDir, LatchFile
from lplots.widgets.workflow import w_workflow

# Ensure spatial coordinates are correctly named
if "X_spatial" not in adata.obsm and "Spatial" in adata.obsm:
    adata.obsm["X_spatial"] = adata.obsm["Spatial"]

# Save locally
local_h5ad = "/tmp/spatial_processed_for_domain_detection.h5ad"
adata.write_h5ad(local_h5ad)

# Upload to Latch
lpath = LPath("latch:///spatial/analysis_output/spatial_processed_for_domain_detection.h5ad")
lpath.upload_from(Path(local_h5ad))

# Configure workflow parameters
params = {
    "input_file": LatchFile("latch:///spatial/analysis_output/spatial_processed_for_domain_detection.h5ad"),
    "run_name": "my_run",
    "output_dir": LatchDir("latch:///Domain_detection_output"),
    "lambda_list": "0.5",
}

# Launch workflow
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
    if res is not None:
        print(res.status)
```

---

## 5. Visualization

Load the output `.h5ad` and visualize domain labels using `w_h5` or Scanpy spatial plots:
- Color by `labels_scaled_gaussian_*`
- Overlay domains on spatial coordinates

---

**Notes**
- Higher lambda → smoother, larger domains  
- Lower lambda → finer, more fragmented domains  
- To compare lambdas, rerun the workflow once per lambda value
