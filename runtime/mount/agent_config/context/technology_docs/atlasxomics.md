# AtlasXomics Analysis Guide

## Analysis Pipeline Overview
1. **Data Loading**
2. **Quality Control** 
3. **Clustering** - Use workflow `wf.__init__.opt_workflow` → `technology_docs/atlasxomics/clustering.md`
4. **Differential Analysis** - Use workflow `wf.__init__.compare_workflow`
5. **Cell Type Annotation** - Use CellGuide marker database → `technology_docs/atlasxomics/cell_type_annotation/overview.md`

## Quality Control

**IMPORTANT**: Skip if clustering requested (workflow includes QC internally).

### Pre-check
Always verify if AnnData is pre-processed: check `adata.obs` for existing QC metrics.

### Key Metrics & Thresholds
Use `snapatac2` library:
```python
import snapatac2 as snap
```

**Adaptive Filtering (per-sample quantiles)**:
- `n_fragments`: [max(q5, 1k), min(q99.5, 50k)]
- `tsse`: ≥ min(q10, 2)
- `frip`: ≥ min(q10, 0.2)
- `nucleosome_signal`: ≤ max(q90, 4)
- `mitochondrial_fraction`: ≤ max(q90, 0.10)

### Metrics

#### 1. Fragment Size Distribution
**Purpose:** Assess nucleosome periodicity and library quality.
**Pattern:** 80-300bp (open chromatin), ~150-200bp (mono), ~300-400bp (di), >500bp (multi/artifacts)
```python
fig = snap.pl.frag_size_distr(data, show=False)
fig.update_yaxes(type="log")
```

#### 2. TSS Enrichment (TSSE)
**Purpose:** Quantify enrichment near transcription start sites.
**Thresholds:** ≥5-10 good, <4 poor
```python
snap.metrics.tsse(data)
```

#### 3. FRiP - Fraction of Reads in Peaks
**Purpose:** Quantify fragments in called peaks (~0.2 good, <0.1 noisy)
```python
snap.metrics.frip(adata, regions, inplace=True, n_jobs=8)
```
**Inputs:**
- `adata`: AnnData or list of AnnData objects
- `regions`: dict mapping peak-set names to BED paths
- `n_jobs`: parallel workers (-1 for all cores)

#### 4. Nucleosome Signal
**Thresholds:** <2 good, >4 over-digested

#### 5. Number of Fragments per Cell
**Field:** `adata.obs["n_fragment"]`
**Thresholds:** <1k dropouts, very high = doublets

#### 6. Mitochondrial Read Fraction
**Field:** `adata.obs["frac_dup"]`
**Thresholds:** >10% = stressed cells, can be relaxed based on adaptive filters for ATAC-seq data.

## Differential Analysis Workflow

`compare_workflow` compares differences in genes, peaks, and motifs between user-defined cluster/condition groupings. Outputs include:

**For genes**:
- `volcano_plot.pdf`: Log2FC vs -log10(p-value) visualization
- `all_genes.csv`: Complete differential expression results from ArchR::getMarkerFeatures
- `marker_genes.csv`: Filtered significant genes with Log2FC, FDR, and enrichment scores

**For peaks**:
- `MA_plot.pdf`: Average signal vs Log2FC visualization
- `all_peaks.csv`: All peak accessibility differences
- `marker_peaks.csv`: Filtered significant peaks with genomic coordinates

**For motifs**:
- `[up/down]Regulated_motifs.csv`: TF motifs ranked by significance
- `[up/down]_enrichment_plot.pdf`: Scatter plots ranked by -log10(FDR)
- `all_motifs.csv`: Complete motif enrichment differences
- `marker_motifs.csv`: Filtered significant motifs

### Parameters
**Required**:
- `project_name` (str)
- `groupings` (LatchFile): JSON with groupA/groupB cell barcodes
- `archrproject` (LatchDir): Path ending in "_ArchRProject"
- `genome` (Enum): "hg38" or "mm10"

### Creating compare_config.json
```python
import json
from pathlib import Path

# Infer comparison column from user request
group_column = "condition"  # or "sample", "cluster"
group_a_value = "Cirrhotic"
group_b_value = "Healthy"

# Get barcodes for each group
group_a_bcs = adata.obs[adata.obs[group_column] == group_a_value].index.tolist()
group_b_bcs = adata.obs[adata.obs[group_column] == group_b_value].index.tolist()

groupings = {
    "groupA": group_a_bcs,
    "groupB": group_b_bcs
}

# Save locally
local_path = "compare_config.json"
with open(local_path, "w") as f:
    json.dump(groupings, f)

# Upload to Latch Data
remote_path = "latch:///compare_config.json"
latch_path = LPath.upload(Path(local_path), remote_path)

# Wrap as LatchFile for workflow input
groupings_file = LatchFile(remote_path)
```

### Launch Comparison
```python
params = {
    "project_name": "my_comparison",
    "groupings": groupings_file,
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

## Data Structure

### Raw Data Paths
**Internal Workspace (13502)**:
- Fragments: `/chromap_outs/[Run_ID]/chromap_output/fragments.tsv.gz`
- Spatial: `/Images_spatial/[Run_ID]/spatial`

**Collaborator Workspaces**:
- Fragments: `.../Raw_Data/[Run_ID]/chromap_output/fragments.tsv.gz`
- Spatial: `.../Raw_Data/[Run_ID]/spatial`

### Workflow Outputs
**Clustering**: `/snap_opts/[project_name]/`
- `setN_*` folders: Each represents different parameter combinations. Contains:
  - `combined.h5ad`: AnnData with .X as tile matrix (genomic bins x cells)
- `figures/`: QC and clustering plots as PDFs to guide parameter selection
- `medians.csv`: QC metrics summary for all samples

**Analysis**: `/snap_outs/[project_name]/` or `.../Processed_Data/[project_name]`
Key files for analysis:
- `combined_sm_ge.h5ad`: **Gene activity scores** for all spots/cells across all samples
  - .X matrix: gene activity (imputed from chromatin accessibility)
  - Used for: gene expression analysis, cell type annotation
- `combined_sm_motifs.h5ad`: **Motif enrichment scores** for all spots/cells
  - .X matrix: TF motif enrichment scores (870 motifs)
  - Used for: transcription factor activity analysis
- `*_ArchRProject/`: ArchR project directory for R-based analysis

Additional outputs:
- `combined_ge.h5ad`: Gene activity without smoothing
- `combined.h5ad`: Original peak/tile matrix
- `[sample]_g_converted.h5ad`: Per-sample gene activity
- `[sample]_m_converted.h5ad`: Per-sample motif enrichment
- `compare_config.json`: Example grouping file for comparisons
- `cluster_coverages/`, `condition_coverages/`, `sample_coverages/`: BigWig coverage tracks
- `figures/`: QC and analysis plots
- `tables/`: Summary statistics

### AnnData Structure
Example: `combined_sm_ge.h5ad` with n_obs × n_vars = 48453 × 24919

**obs** (cell/spot metadata):
- `n_fragment`: Total fragments per cell
- `frac_dup`: Duplication rate
- `frac_mito`: Mitochondrial fraction
- `tsse`: TSS enrichment score
- `log10_frags`: Log10 of fragment count
- `on_off`: Spot active (1) or empty (0)
- `row`, `col`: Array grid coordinates
- `xcor`, `ycor`: Physical spatial coordinates
- `sample`: Sample/run identifier (e.g., "D02297")
- `condition`: Experimental condition
- `cluster`: Cluster assignment from analysis

**obsm** (embeddings):
- `X_umap`: UMAP coordinates for visualization
- `spatial`: Spatial coordinates array

**obsp** (graphs):
- `spatial_connectivities`: Spatial neighbor graph
- `spatial_distances`: Distance matrix between neighbors
