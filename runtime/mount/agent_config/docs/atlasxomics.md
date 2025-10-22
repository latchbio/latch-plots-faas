## Analysis Guideline

This is the **authoritative step-by-step pipeline** for AtlasxOmics experiment. Follow steps in order. 

1. **Experiment Setup** - If not clear from original request, ask users to confirm if they want to perform analysis on **gene activity score AnnData** (recommended) or **motif enrichment scores AnnData**. 
2. **Data Loading** - load data using **Scanpy** and display it with `w_h5`.
3. **Clustering (workflow only)** - Launch the AtlasXOmics clustering workflow using `w_workflow(wf_name="wf.__init__.opt_workflow", ...)`. Fallback to `scanpy` only if this fails.
4. **Differential Gene Activity or Motif Enrichment Comparison** - Use `w_workflow(wf_name="wf.__init__.compare_workflow", ...)`
5.  **Cell Type Annotation** - assign biological meaning to clusters using gene sets. 

The section below defines detailed guidelines for each of the above steps.

### **Workflow Rules**

<workflow_rules>
**HARD RULES — DO NOT VIOLATE:**
1. **Forms only**. Collect every workflow input via `lplots.widgets`. Do not ask questions in prose.
2. **No hardcoded params**. Every workflow param must come from a widget value or a derived default shown in a widget.
3. **Workflow input paths are remote**. Never pass local paths; upload then use `latch://` in `LatchFile`/`LatchDir`.
4. **Exact workflow version**. Use the exact string provided in examples; do not change or omit it.
5. Must end with `execution = w.value` so a button is displayed to run the workflow.
</workflow_rules>

### **Data Loading**: 
- Identify the correct node ID for either `combined_sm_ge.h5ad` (gene activity scores) or `combined_sm_motifs.h5ad` (motif enrichment), and load the file using `LPath("latch://<node_id>.node")`
- After loading, always display the AnnData object with the `w_h5` widget.

### **Clustering (workflow only)**: 
- Use `w_workflow` to launch the AtlasXOmics clustering workflow.
- If the user want to cluster only a subset of the AnnData, include the optional parameter `adata_subsetted_file` (as a LatchFile) in params. Otherwise, omit it entirely.

### Clustering Workflow Parameters

- Strictly follow <workflow_rules>.

#### **Required**
- `project_name` *(str)*  
  Name for the output folder.

- `genome` *(Enum)*  
  One of: `"hg38"`, `"mm10"`, `"rnor6"`.

- `runs` *(List[Run])*  
  List of individual sample inputs, each with:
  - `run_id` *(str)* — Sample identifier  
  - `fragments_file` *(LatchFile)* — BED-like file of ATAC-seq fragments  
  - `spatial_dir` *(LatchDir)* — Tissue image + metadata folder  
  - `condition` *(str, optional)* — e.g., `"control"` or `"diseased"`

#### **Optional**
- `adata_subsetted_file` *(LatchFile)*  
  Subset AnnData to cluster (omit if clustering full object)

- `tile_size` *(int)*  
  Genomic bin size (default: `5000`)

- `n_features` *(List[int])*  
  Top accessible tiles to use, e.g., `[25000]`

- `resolution` *(List[float])*  
  Clustering resolution, e.g., `[1.0]`

- `varfeat_iters` *(List[int])*  
  Iterations for variable feature selection, e.g., `[1]`

- `n_comps` *(List[int])*  
  Dimensionality reduction components, e.g., `[30]`

- `min_cluster_size` *(int)*  
  Minimum cells per cluster

- `min_tss` *(float)*  
  Minimum TSS enrichment score per cell

- `min_frags` *(int)*  
  Minimum fragments per cell

- `pt_size` *(int or None)*  
  Point size override for spatial plots

- `qc_pt_size` *(int or None)*  
  Point size override for QC spatial plots

Use the code below as a template. The code will generate a "Launch" button.
Instruct users to click this button to start the workflow.
                     
#### Example Implementation

```python
from dataclasses import dataclass
from latch.types import LatchFile, LatchDir
from latch.ldata.path import LPath
from lplots.widgets.ldata import w_ldata_picker
from lplots.widgets.text import w_text_output, w_text_input

class Genome(Enum):
    mm10 = "mm10"
    hg38 = "hg38"
    rnor6 = "rnor6"

@dataclass
class Run:
    run_id: str
    fragments_file: LatchFile
    spatial_dir: LatchDir
    condition: str = "None"

# 1) Expose widget for user to pick the top-level remote folder (returns LPath)
raw_data_dir = w_ldata_picker(label="Raw Data Directory")
if raw_data_dir.value is None:
    w_text_output(
        content="Select the top-level folder containing sample subfolders.",appearance={"message_box": "warning"}
    )
    exit(0)

root: LPath = raw_data_dir.value  # already an LPath — use directly

# 2) Samples + condition map
samples = adata.obs["sample"].unique()
cond_map = (adata.obs.groupby("sample")["condition"].first().to_dict()
            if "condition" in adata.obs.columns
            else {s: "None" for s in samples})

# 3) Helper to join remote paths (LPath uses "/" join)
def rpath(*parts: str) -> LPath:
    p = root
    for part in parts:
        p = p / part
    return p  # still LPath

# 4) Build runs (convert to LatchFile/Dir via .path, not str())
runs = []
for s in samples:
    print(s)
    frag_lp: LPath = rpath(str(s), "chromap_output", "fragments.tsv.gz")
    spat_lp: LPath = rpath(str(s), "spatial")

    runs.append(
        Run(
            run_id=str(s),
            fragments_file=LatchFile(frag_lp.path),  # latch://...
            spatial_dir=LatchDir(spat_lp.path),      # latch://...
            condition=cond_map[str(s)],
        )
    )

# 5) Construct params with minimal form for user to input project name, genome, resolution

project_name = w_text_input(label="Project name", default="atlasx_clustering")
genome = w_select(label="Genome", options=[g.value for g in Genome], default="hg38")
resolution = w_text_input(label="Resolution(s)", default="0.5")

def to_list_of_floats(text: str):
    return [float(x.strip()) for x in text.split(",") if x.strip()]

params = {
     "runs": runs, 
     "adata_subsetted_file": LatchFile("latch:///adata_subset.h5ad"),
     "genome": genome.value,
     "project_name": project_name.value,
     "tile_size": 5000,
     "n_features": [25000],
     "resolution": to_list_of_floats(resolution.value),
     "varfeat_iters": [1],
     "n_comps": [30],
     "min_cluster_size": 20,
     "min_tss": 2.0,
     "min_frags": 10,
     "pt_size": None,
     "qc_pt_size": None,
}

w = w_workflow(
  wf_name="wf.__init__.opt_workflow",
  version="0.3.5-9e16e4",
  params=params,
  label="Run clustering workflow",
)

execution = w.value
```

#### Example Implementation
```python
viewer = w_h5(ann_data=adata)
value = viewer.value

# Get the current selections
if value['lasso_points']:
    print(f"User selected {len(value['lasso_points'])} regions")
    print(f"The embedding used for the lasso selection is {value['lasso_points_obsm']}")
    for i, region in enu merate(value['lasso_points']):
        print(f"Region {i}: {len(region)} points")

# Proceed to create an `adata_subset` based on lasso-selected points
```

### **Differential Gene Activity or Motif Enrichment Comparison (workflow only)**
- Use `w_workflow` to launch the AtlasXOmics comparison workflow.  
- Automatically infer the correct grouping column from `adata.obs` (`condition`, `sample`, or `cluster`).  
- Programmatically generate and upload a `compare_config.json` file as a `LatchFile` for workflow input.
- Strictly follow <workflow_rules>.

### Comparison Workflow Parameters

#### **Required**
- `project_name` *(str)*  
  Name for the output folder.

- `groupings` *(LatchFile)*  
  A JSON file (`compare_config.json`) containing two keys:  
  - `groupA`: list of AnnData cell barcodes for the first group  
  - `groupB`: list of AnnData cell barcodes for the second group  

- `archrproject` *(LatchDir)*  
  Path to an ArchR project directory on Latch Data (usually ends with `_ArchRProject`).

- `genome` *(Enum)*  
  Either `"hg38"` or `"mm10"`.

#### Example Implementation

```python
from latch.types import LatchFile, LatchDir
from lplots.widgets.workflow import w_workflow

params = {
    "project_name": "my_comparison",
    "groupings": LatchFile("latch:///compare_config.json"),
    "archrproject": LatchDir("latch:///Kostallari_SOW313_ATAC_ArchRProject"),
    "genome": "hg38",
}

w = w_workflow(
    wf_name="wf.__init__.compare_workflow",
    version="0.7.1-8484d6-wip-4ae938",
    params=params,
    label="Launch Comparison Workflow"
)

execution = w.value
```

#### How to construct `compare_config.json`
To run a comparison workflow, you must generate a `compare_config.json` file that defines which cells belong to each group.

What to Ask the User:
- Comparison Target: Ask the user what biological groups or conditions they want to compare. For example:
    - "Diseased vs Healthy"
    - "Sample A vs Sample B"
    - "Cluster 5 vs Cluster 7"

How to Identify Groups:
- You must infer which column in adata.obs encodes this grouping (e.g. "condition", "sample", or "cluster").
- Then, for each group (A and B), filter all cell barcodes in `adata.obs.index` that match the selected value.

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

### **Cell Type Annotation**

- If the dataset context is unclear, first ask the user to confirm the **organism** and **tissue type**. **Do NOT proceed** until the user has answered the question. 
- **Always render a form with sensible defaults** to avoid tedious manual input. The form should support **multiple candidate cell types**, e.g. one row per cell type:
  - `cell_type`: **text input**, pre-filled with a common or inferred cell type.
  - `marker_genes`: **multiselect widget**, pre-populated with default marker genes for that cell type.
- You **must auto-populate all fields with reasonable defaults using domain knowledge**. Users should only adjust values if needed, not enter them from scratch.
- Add a **button** after the form to trigger gene set scoring. 
- For **spatial ATAC-seq** only, infer cell identity by computing **gene activity or gene set scores** (e.g., `scanpy.tl.score_genes`) and ranking cell types based on marker enrichment.

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
- Folder ending with `_coverages`: Contains `.bw` files which can be visualized in the IGV browser.

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
