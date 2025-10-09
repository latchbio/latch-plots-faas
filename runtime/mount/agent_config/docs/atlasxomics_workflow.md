## Analysis Guideline

This is the **authoritative step-by-step pipeline** for AtlasxOmics experiment. Follow steps in order. 

1. **Data Loading** - load data using **Scanpy** and `w_h5`.
2. **Experiment Setup** - Ask users to confirm if they want to perform analysis on **gene activity score AnnData** (recommended) or **motif enrichment scores AnnData**. 
3. **Clustering (workflow only)** - Confirm if clustering applies to all cells or a subset (subset = a single sample, condition, or lasso-selected region). Then launch the AtlasXOmics clustering workflow using `w_workflow(wf_name="wf.__init__.opt_workflow", ...)` **Never use scanpy.tl.leiden** or any ad-hoc clustering. 
9. **Differential Gene Activity or Motif Enrichment Comparison** - Use `w_workflow(wf_name="wf.__init__.compare_workflow", ...)`
10. **Cell Type Annotation** - assign biological meaning to clusters using gene sets. 

The section below defines detailed guidelines for each of the above steps. 

### **Data Loading**: 
- Load the motif and gene H5AD into two adata objects and view them with `w_h5` widget. 
- The structure of the AtlasXOmics input folder is standard, so assume and check that `combined_sm_ge.h5ad` and `combined_sm_motifs.h5ad` exist, without asking users to confirm. 

### **Clustering (workflow only)**: 
- If users want to perform clustering, launch the AtlasXOmics clustering workflow using the `w_workflow` widget. 
- The workflow requires precise user input because each field maps directly to workflow parameters in the code. **Always output a form with widgets for users to fill out required inputs**. 
- You must **parse user answers**, **normalize them into the required formats**, and then construct the `params` dictionary exactly as shown in the example.  
- If the user specifies they want to cluster only a subset of the AnnData, include the optional parameter `adata_subset` (as a LatchFile) in params. Otherwise, omit it entirely.

#### Required User Inputs
- `run_id`: An identifier for the Run.
- `fragments_file`: A file on Latch Data. A BED-like, tab-delimited file in which each row contains an ATAC-seq fragment.
- `condition` (optional): An experimental Condition descriptor (i.e., “control”, 'diseased’).
- `spatial_dir`: A directory on Latch Data, containing tissue images and experiment metadata.
- `genome`: Can either be "hg38" for humans or "mm10" for mouse 
- `adata_subset (optional)`: A LatchFile representing the subset of the AnnData to cluster.

Use the code below as a template. The code will generate a “Launch” button.
In your summary response, explicitly instruct users to click this button to start the workflow.
                     
#### Example Implementation

```python
from lplots.widgets.workflow import w_workflow
from latch.types import LatchFile, LatchDir
from lplots.widgets.workflow import w_workflow
import sys
import importlib.util

# Import statements
utils_path = "/opt/latch/plots-faas/optimize_snap/wf/utils.py"
spec = importlib.util.spec_from_file_location("wf_utils", utils_path)
wf_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wf_utils)

Genome = wf_utils.Genome
Run = wf_utils.Run

params = {
    "project_name": "Test Project",
    "runs": [
        Run(
            run_id="D02301_NG07297",
            fragments_file=LatchFile("latch://38438.account/chromap_output/fragments.tsv.gz"),
            condition="Healthy",
            spatial_dir=LatchDir(
                "latch://38438.account/spatial"
            )
          )
    ],
    "adata_subset": LatchFile("latch://38438.account/adata_subset.h5ad"), 
    "genome": Genome["hg38"],
    "resolution": [1.0],
    "tile_size": 5000, # default
    "n_features": [25000], # default
    "varfeat_iters": [1], # default
    "n_comps": [30], # default
    "min_cluster_size": 20, # default
    "min_tss": 2.0, # default
    "min_frags": 10 # default
}

w = w_workflow(
    wf_name="wf.__init__.opt_workflow",
    version="0.3.5-f10e3c-wip-9a1d7b",
    params=params,
    label="Run Clustering Workflow",
)
execution = w.value
```

### **Lasso Selection Tasks or AnnData Subset Selection Tasks**
- Use when users want to perform clustering on specific spatial regions they've selected
- Reference Latch Plots documentation to retrieve lasso selected values and perform targeted clustering

#### Example Implementation
```python
viewer = w_h5(ann_data=adata)
value = viewer.value
if value.lasso_points:
    print(f"User selected {len(value.lasso_points)} regions")
    print(f"The embedding used for the lasso selection is {value.lasso_points_obsm}")
# Proceed to create an `adata_subset` based on lasso-selected points
```

### Differential Gene Activity or Motif Enrichment Comparison Workflow

#### Required User Inputs: 
**Always create a form with widgets to collect user inputs**:
- `project_name`: a string for project name
- `groupings`: A `LatchFile` named `compare_config.json` that stores a dictionary with two keys: `groupA` and `groupB`. Each key has a list of AnnData cell barcodes as values. 
- `archrprokect`: A `LatchDir` that stores an ArchR project on Latch Data. The folder often ends with `_ArchRProject`. 

```python
params = {
        "project_name": "USER_INPUT",
        "groupings": LatchFile(remote_bcs.path),
        "archrproject": LatchDir(archrproj_dir.path),
        "genome": "hg38" # Can be "hg38" or "mm10"
    }
    
w = w_workflow(
    wf_name="wf.__init__.compare_workflow",
    version="0.7.1-8484d6-wip-4ae938",
    params=params,
    label="Launch Workflow"
)

execution = w.value
```

**Constructing groupings for Differential Comparison**
To run a comparison workflow, you must generate a `compare_config.json` file that defines which cells belong to each group.

What to Ask the User:
- Comparison Target: Ask the user what biological groups or conditions they want to compare. For example:
    - "Diseased vs Healthy"
    - "Sample A vs Sample B"
    - "Cluster 5 vs Cluster 7"

How to Identify Groups:
- You must infer which column in adata.obs encodes this grouping (e.g. "condition", "sample", or "cluster").
- Then, for each group (A and B), filter all cell barcodes in adata.obs.index that match the selected value.

```python
import json
from pathlib import Path

group_column = "condition"  # inferred from context or user input
group_a_value = "Cirrhotic"
group_b_value = "Healthy"

group_a_bcs = adata.obs[adata.obs[group_column] == group_a_value].index.tolist()
group_b_bcs = adata.obs[adata.obs[group_column] == group_b_value].index.tolist()

groupings = {
    "groupA": group_a_bcs,
    "groupB": group_b_bcs
}

local_path = "compare_config.json"

with open(local_path, "w") as f:
    json.dump(groupings, f)

# Upload to Latch Data
remote_path = "latch:///compare_config.json"
latch_path = LPath.upload(Path(local_path), remote_path)

# Wrap as LatchFile for workflow input
groupings_file = LatchFile(remote_path)
```

## Data Assumptions

AtlasXOmics datasets follow a standard output folder schema across all projects.
Assume this structure exists without asking the user to confirm.

### Example Folder Layout
root/
├── cluster_coverages/
├── condition_coverages/
├── figures/
├── Kostallari_SOW313_ATAC_ArchRProject/
├── sample_coverages/
├── tables/
│
├── combined_ge.h5ad
├── combined.h5ad
├── combined_motifs.h5ad
├── combined_sm_ge.h5ad
├── combined_sm_motifs.h5ad
│
├── compare_config.json
│
├── D02297_NG07294_g_converted.h5ad
├── D02297_NG07294_m_converted.h5ad
├── D02297_NG07294_SeuratObjMotif.rds
├── D02297_NG07294_SeuratObj.rds
│
├── D02310_NG07345_g_converted.h5ad
├── D02310_NG07345_m_converted.h5ad
├── D02310_NG07345_SeuratObjMotif.rds
├── D02310_NG07345_SeuratObj.rds
│
├── enrichMotifs_clusters.rds
├── enrichMotifs_condition_1.rds
├── enrichMotifs_sample.rds
│
├── markersGS_clusters.rds
├── markersGS_condition_1.rds
├── markersGS_sample.rds


**Important files to pay attention to**
- `combined_sm_ge.h5ad`: Stores gene activity score for every spot
- `combined_sm_motifs.h5ad`: Stores motif enrichment score for every spot
- Both files follow the standard AnnData structure with `obs`, `uns`, `obsm`, and `obsp` components as detailed below. 
- Folder ending with `_ArchRProject`: An `ArchR` project 

**Standard Fields inside an AnnData Object**
Example structure of `combined_sm_ge.h5ad`
AnnData object with n_obs × n_vars = 48453 × 24919
    obs: 'n_fragment', 'frac_dup', 'frac_mito', 'on_off', 'row', 'col', 'xcor', 'ycor', 'sample', 'condition', 'tsse', 'log10_frags', 'cluster'
    obsm: 'X_umap', 'spatial'
    obsp: 'spatial_connectivities', 'spatial_distances'

Example structure of `combined_sm_motifs.h5ad`
AnnData object with n_obs × n_vars = 48453 × 870
    obs: 'n_fragment', 'frac_dup', 'frac_mito', 'on_off', 'row', 'col', 'xcor', 'ycor', 'sample', 'condition', 'tsse', 'log10_frags', 'cluster'
    obsm: 'X_umap', 'spatial'

- **Always** assume these fields exist unless otherwise specified.
- Never ask the user to confirm presence of fields such as `condition`, `sample`, or `cluster` — assume they exist in `adata.obs`.
