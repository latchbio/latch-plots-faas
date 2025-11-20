# Clustering Workflow

`optimize_snap` is the LatchBio Workflow for assessing DBiT-seq epigenomic experiment quality and systematically exploring clustering parameter combinations. It generates multiple H5AD outputs—each clustered with different parameter settings—along with summary statistics to guide downstream analysis.

**ABSOLUTE RULE**: The Workflow contains **rigorously tested, gold-standard SnapATAC2 code internally and should ALWAYS be ysed OVER writing code from scratch.**

---

### Input Data & Decision Logic

**The clustering workflow ALWAYS requires raw data folders (for the `runs` parameter).**

**User scenarios:**

1. **Provides raw data only** → Standard clustering
   - Validate: Each sample needs `chromap_output/fragments.tsv.gz` + `spatial/`
   - Launch workflow with validated samples

2. **Has AnnData subset** → Subset clustering  
   - Save: `adata_subset.to_h5ad("subset.h5ad")`
   - Ask: "Where are the raw data folders?"
   - Launch with subset file + raw data paths

3. **Has clustered H5AD** → Clarify intent
   - Ask: "Keep current clustering or re-cluster?"
   - If re-cluster → Ask: "Where are the raw data folders?"
   - Launch with raw data paths

**Key: Raw data is mandatory - the workflow cannot run without fragment files and spatial folders.**
---

### Parameters
**Required**:
- `project_name` (str)
- `genome` (Enum): "hg38", "mm10", "rnor6"
- `runs` (List[Run]): Each with:
  - `run_id` (str)
  - `fragments_file` (LatchFile)
  - `spatial_dir` (LatchDir)
  - `condition` (str, optional)

**Recommended Defaults**:
- `n_features`: [25_000, 50_000, 100_000]
- `resolution`: [0.5, 1.0, 1.25, 1.5] (AtlasXomics needs ≥0.5)
- `n_comps`: [30, 50]

**Optional Parameters**:
- `adata_subsetted_file` (LatchFile): Subset AnnData to cluster
- `tile_size` (int): Genomic bin size (default: 5000)
- `varfeat_iters` (List[int]): Variable feature iterations (default: [1])
- `min_cluster_size` (int): Minimum cells per cluster
- `min_tss` (float): Minimum TSS enrichment
- `min_frags` (int): Minimum fragments per cell
- `pt_size` (int or None): Point size for spatial plots
- `qc_pt_size` (int or None): Point size for QC plots

### Implementation Template for Collecting User Inputs
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
        content="Select the top-level folder containing sample subfolders.",
        appearance={"message_box": "warning"}
    )
    exit(0)

root: LPath = raw_data_dir.value  # already an LPath

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

# 5) Construct params with minimal form for user input
project_name = w_text_input(label="Project name", default="atlasx_clustering")
genome = w_select(label="Genome", options=[g.value for g in Genome], default="hg38")
resolution = w_text_input(label="Resolution(s)", default="0.5")

def to_list_of_floats(text: str):
    return [float(x.strip()) for x in text.split(",") if x.strip()]
```

### Implementation Template for Launching the Workflow and Waiting for it to Complete

```python
params = {
     "runs": runs, 
     "adata_subsetted_file": LatchFile("latch:///adata_subset.h5ad"),  # optional
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
  key="clustering_workflow_run_1",
  version="0.3.5-9e16e4",
  params=params,
  automatic=True,
  label="Run clustering workflow",
)

execution = w.value

if execution is not None:
  res = await execution.wait()

  if res is not None and res.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
      # inspect workflow outputs for downstream analysis
      workflow_outputs = list(res.output.values())
```

---

### Self-evaluation Checklist

Before launching the workflow, evaluate your work:

- Create widgets: w_ldata_picker, w_text_input, w_select with defaults shown
- Verify all LatchFile/LatchDir paths exist on Latch
- Resolution ≥ 0.5
- For every sample, validate the following file/ folder exists: 
    - Fragment files: /[sample]/chromap_output/fragments.tsv.gz
    - Spatial dirs: /[sample]/spatial/
- If there is AnnData object in notebook, confirm samples match between adata and raw data
