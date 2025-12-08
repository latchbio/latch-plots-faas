# AtlasXomics Analysis Guide

<atlasx_analysis_overview>

## Analysis Pipeline Overview
1. **Data Loading**
2. **Quality Control** → `technology_docs/atlasxomics/qc.md`
3. **Clustering** - Use workflow `wf.__init__.opt_workflow` → `technology_docs/atlasxomics/clustering.md`
4. **Differential Analysis** - Use workflow `wf.__init__.compare_workflow` → `technology_docs/atlasxomics/de.md` 
5. **Cell Type Annotation** - Use CellGuide marker database → `technology_docs/atlasxomics/cell_type_annotation/overview.md`

</atlasx_analysis_overview>

---

<pre_analysis_step>

MANDATORY: Invoke the `redeem_package` tool to install required AtlasXomics tools into the workspace.
  - `package_code`: `2428814b149447a4c354b3cb4520095b77955bf99cb3eedfef20b920a2a7d3d7`
  - `package_version_id`: `405`

</pre_analysis_step>

</atlasx_data_structure>

## Data Structure

### Raw Data Paths
**Internal Workspace (13502)**:
- Fragments: `/chromap_outs/[Run_ID]/chromap_output/fragments.tsv.gz`
- Spatial: `/Images_spatial/[Run_ID]/spatial`

**Collaborator Workspaces**:
- Fragments: `.../Raw_Data/[Run_ID]/chromap_output/fragments.tsv.gz`
- Spatial: `.../Raw_Data/[Run_ID]/spatial`

### Workflow Outputs

All workflow outputs are stored under the root Latch Data directory (`latch://<account_id>.id/` or `latch:///`):

**Differential Expresion** `/compare_outs/[project_name]/`
- **`gene_results/`**
  - Volcano plot  
  - `all_genes.csv` — all gene-level test results  
  - `marker_genes.csv` — significant genes after filtering  

- **`peak_results/`**
  - MA plot  
  - `all_peaks.csv` — all peak-level test results  
  - `marker_peaks.csv` — significant peaks after filtering  

- **`motif_results/`**
  - `[up/down]Regulated_motifs.csv` — motifs ranked by significance  
  - `[up/down]_enrichment_plot` — motif enrichment scatter plot (`-log10(FDR)`)  
  - `all_motifs.csv` — all motif-level test results  
  - `marker_motifs.csv` — significant motifs after filtering  

- **`coverages/`**
  - BigWig coverage tracks for visualization  

- **`[project_name]_ArchRProject/`**
  - Updated ArchRProject with differential results attached  


**Clustering**: `/snap_opts/[project_name]/`
- `setN_*` folders: Each represents different parameter combinations. Contains:
  - `combined.h5ad`: AnnData with .X as tile matrix (genomic bins x cells)
- `figures/`: QC and clustering plots as PDFs to guide parameter selection
- `medians.csv`: QC metrics summary for all samples

**Analysis**: `/snap_outs/[project_name]/` or `.../Processed_Data/[project_name]`
Key files for analysis:
- `combined_sm_ge.h5ad`: **Gene activity scores** for all spots/cells across all samples. **Recommended for analyses**
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

</atlasx_data_structure>
