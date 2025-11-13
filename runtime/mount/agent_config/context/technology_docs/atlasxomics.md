## Analysis Guideline

1. **Data Loading**
2. **Quality Control**
3. **Clustering (workflow only)** - Launch the AtlasXomics workflow for quality control and clustering using `w_workflow(wf_name="wf.__init__.opt_workflow", ...)`. This contains **gold-standard, well-tested code** that uses SnapATAC2 internally, and should _always_ be used over writing code from scratch. 
4. **Differential Gene Activity or Motif Enrichment Comparison** - Use `w_workflow(wf_name="wf.__init__.compare_workflow", ...)`
5. **Cell Type Annotation** — Use CellGuide marker database (see file `technology_docs/marker_cell_typing.md`)

The section below defines detailed guidelines for each of the above steps.

### **Quality Control**:

**IMPORTANT**: Skip this step if user requests clustering because the clustering workflow contains internal code for QC and filtering. 

- Use the `snapatac2` library for computing and visualizing ATAC-seq quality metrics.
- **MANDATORY**: Before running any QC or filtering steps, **verify whether the AnnData object has already been pre-processed or quality-controlled.** If yes, skip QC and ask users if they actually want to re-perform QC. 

- Use `snapatac2` library for computing and visualizing all ATAC-seq QC metrics.
- Always check first whether QC metrics already exist in `adata.obs`. 
- **If QC already exists and passes reasonable thresholds**: skip recomputation and filtering; optionally confirm with the user if they want to re-run QC and filtering. 
- **If QC is missing:** compute required metrics (fragment size distribution, TSSE, FRiP, nucleosome signal, fragments per cell, mitochondrial read fraction), and proceed to filtering. 

```python
import snapatac2 as snap
```
- **Key QC metrics**: Check which metrics already exist in `adata.obs`, run adaptive filtering with those first, and only compute any missing metrics afterward if needed. 
  - Fragment Size Distribution
  - TSS Enrichment (TSSE)
  - FRiP — Fraction of Reads in Peaks
  - Nucleosome Signal
  - Number of Fragments per Cell
  - Mitochondrial Read Fraction

- **Adaptive, per-sample QC filtering**:
  Inputs (adata.obs): n_fragments, tsse, frip, nucleosome_signal, mitochondrial_fraction
  Batch key: sample
  Heuristic (per-batch quantiles):
    n_fragments: keep [max(q5, 1k), min(q99.5, 50k)]
    tsse: ≥ min(q10, 2)
    frip: ≥ min(q10, 0.2)
    nucleosome_signal: ≤ max(q90, 4)
    mitochondrial_fraction: ≤ max(q90, 0.10)

#### 1. Fragment Size Distribution

**Purpose:** Assess nucleosome periodicity and library quality.  
**Expected pattern:**
- **80–300 bp:** Nucleosome-free (open chromatin)
- **~150–200 bp:** Mono-nucleosome peak
- **~300–400 bp:** Di-nucleosome peak
- **>500 bp:** Multi-nucleosome or artifacts

**Example:**
```python
fig = snap.pl.frag_size_distr(data, show=False)
fig.update_yaxes(type="log")
```

#### 2. TSS Enrichment (TSSE)

**Purpose:** Quantify enrichment of accessible fragments near transcription start sites.  
- **High TSSE (≥ 5–10):** Strong promoter accessibility, good quality  
- **Low TSSE (< 4):** Poor signal, low complexity, or over-digestion

**Example:**
```python
snap.metrics.tsse(data)
```

#### 3. FRiP — Fraction of Reads in Peaks

**Purpose:** Quantify the share of fragments falling inside called peaks; higher FRiP means cleaner regulatory signal (≈0.2 good, <0.1 noisy).

**Example:**
```python
snap.metrics.frip(adata, regions, inplace=True, n_jobs=8)
```

**Inputs:**  
- `adata`: AnnData or list of AnnData objects to annotate; writes scores to `adata.obs` when `inplace=True`.  
- `regions`: dict mapping peak-set names to BED paths or genomic interval lists.  
- `n_jobs`: parallel workers (use `-1` for all cores).

**Note:** Run `snap.pp.import_data(...)` beforehand to load fragment data.

#### 4. Nucleosome Signal

**Purpose:** Ratio of mono/di-nucleosomal to short fragments.  
- **Low (< 2):** Good chromatin accessibility  
- **High (> 4):** Over-digested or low-quality libraries

#### 5. Number of Fragments per Cell (`adata.obs["n_fragment"]`)

**Purpose:** Assess sequencing depth and data sparsity per cell/barcode.  
- **Low fragments (< 1 k):** Dropouts or ambient noise  
- **Extremely high:** Doublets or multiplets

#### 6. Mitochondrial Read Fraction (`adata.obs["frac_dup"]`)

**Purpose:** Detect low-quality or dying cells with excessive mitochondrial reads.  
- **High (> 10 %):** Possible cell stress or broken nuclei

---

### **Clustering (workflow only)**: 

`optimize_snap` is the LatchBio Workflow for assessing DBiT-seq epigenomic experiment quality and systematically exploring clustering parameter combinations. Given fragment files from a single-cell ATAC-seq preprocessing pipeline (e.g., Chromap) and spatial coordinates, it generates multiple H5AD outputs—each clustered with different parameter settings—along with summary statistics to guide downstream analysis.

This Workflow provides a Python-based alternative to the R-based ArchR pipelines, using SnapATAC2 and Scanpy. As reported in recent literature and verified internally (e.g., Luo 2024), SnapATAC2 can produce clustering results that differ from ArchR, sometimes improving biological resolution.

**IMPORTANT**: The Workflow contains **rigorously tested, gold-standard SnapATAC2 code and should always be preferred over writing code from scratch in the notebook.** 

- Launch the AtlasXomics clustering workflow using `w_workflow`.
- To cluster only a subset of the AnnData, supply the optional `adata_subsetted_file` parameter (as a LatchFile). Otherwise, leave it out entirely.

### Clustering Workflow Parameters

- Because clustering outcomes are sensitive to input settings, always set default `n_features`, `resolution`, and `n_comps` to multiple values so users can compare results and select the most meaningful clustering after the workflow finishes running. 
- **Recommended defaults for AtlasXomics data:**
  - `n_features`: [25_000, 50_000, 100_000]
  - `resolution`: [0.5, 1.0, 1.25, 1.5] — AtlasXomics datasets typically need ≥ 0.5 even for coarse clustering; lower values almost always collapse to a single cluster.
  - `n_comps`: [30, 50]

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
  Top accessible tiles to use, e.g., `[25000, 50000, 100000]`

- `resolution` *(List[float])*  
  Clustering resolution, e.g., `[0.5, 1.0, 1.25, 1.5]`

- `varfeat_iters` *(List[int])*  
  Iterations for variable feature selection, e.g., `[1]`

- `n_comps` *(List[int])*  
  Dimensionality reduction components, e.g., `[30, 50]`

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

### **Differential Gene Activity or Motif Enrichment Comparison (workflow only)**

compare clusters is a latch.bio workflow for comparing differences in genes, peaks, and motifs between user-defined cluster/condition groupings within an ArchRProject. Provided an ArchRProject and grouping specifications, compare clusters generates,

- for genes
  - volcano plot
  - all_genes.csv: all genes with test results from ArchR::getMarkerFeatures
  - marker_genes.csv: all_genes filtered and scored with significance thresholds; data for the volcano plot.

- for peaks
  - MA Plot
  - all_peaks: all peaks with test results from ArchR::getMarkerFeatures
  - marker_peaks.csv: all_genes filtered and scored with significance thresholds; data for the volcano plot.

- for motifs
  - [up/down]Regulated_motifs.csv: Up or down regulated motifs ranked by significance values; see ArchR docs
  - [up/down] enrichment plot: scatter plot of motifs ranked by -log10(FDR)
  - all_motifs.csv: all motifs with test results from ArchR::getMarkerFeatures
  - marker_motifs.csv: all_motifs filtered and scored with significance thresholds; data for the volcano plot.

- Use `w_workflow` to launch the AtlasXomics comparison workflow.  
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

- Comparison Target: Ask the user what biological groups or conditions they want to compare.
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

## Data Assumptions

AtlasXomics datasets have the following data conventions. Assume this structure exists without asking the user to confirm.

### Raw Data
Raw data consists of fragment files (fragments.tsv.gz) and 'spatial/' folders which contain images, image metadata, and barcode-image mappings stored as csv files. Every experiment (designated with the unique 'Run ID' Dxxxxx where x is a digit), is associated with distinct raw data. Fragment files and spatial folders are stored on different files paths in Latch Data.

In the default AtlasXomics Workspace (13502), fragment files are stored in the path `/chromap_outs/[Run_ID]/chromap_output/fragments.tsv.gz`.  Spatial folders are stored in the path `/Images_spatial/[Run_ID]/spatial`. 

In collaborator Workspaces (not 13502), fragment files and spatial folders are stored together in the parent directory corresponding to Run ID.  Frament files are store in the path `.../Raw_Data/[Run_ID]/chromap_output/fragments.tsv.gz`, BED file at `.../Raw_Data/[Run_ID]/chromap_output/aln.bed`, spatial folders at `.../Raw_Data/[Run_ID]]/spatial`

### Workflow Outputs

#### Clustering Workflow (optimize_snap)
The clustering workflow (optimize_snap, wf.__init__.opt_workflow) stores outputs on the path `/snap_opts/[project name]/` where "project name" is designated by the user in the Workflow input parameters. The output directory has the following structure:
[project name]/
├── figures/
├── medians.csv
├── set1_ts5000-vf500000-cr1-0-vi1-nc30
├── set2_ts5000-vf500000-cr1-0-vi1-nc40

Each folder with the prefix 'setN' corresponds to a combination of input clustering parameters.  Each contains a 'combined.h5ad' file which stores the AnnData object generated with the specified clustering parameters, with .X as a tile matrix.

The figures/ directory contains plots saved as .pdf files for QC and to guide selection of cluster parameters for downstream analysis. The medians.csv contains QC metrics for the project.

#### atx_snap and create ArchRProject Workflows

These two Workflows create files to be analyzed in Plots.  A Workflow takes as input one or multiple Run IDs and corresponding raw data.  It creates the outputs detailed below.  The `combined_sm_ge.h5ad` and `combined_sm_motifs.h5ad` files are AnnData objects with .X as gene accessibility data and motif enrichment data, respectively.  Each object contains data for all Run IDs specified in the inputs. Runs are designated by the column 'sample' in AnnData .obs.

In the default AtlasXomics Workspace (13502), results are stored in `/snap_outs/[project name]` where project name is specified in the Execution inputs.  In customer Workspaces, data is stored in `.../Processed_Data/[project name]`.

project_name/
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

**Important files to pay attention to**
- `combined_sm_ge.h5ad`: Stores gene activity score for every spot for all Run IDs
- `combined_sm_motifs.h5ad`: Stores motif enrichment score for every spot for all Run IDs
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
