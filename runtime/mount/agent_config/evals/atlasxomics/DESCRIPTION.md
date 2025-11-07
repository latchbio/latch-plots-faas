# AtlasXomics Spatial ATAC-seq Workflow + Evals

The sections below summarise the evaluations that live **inside**:

- `atlasxomics/qc`
- `atlasxomics/clustering`
- `atlasxomics/cell_typing`
- `atlasxomics/de`

Each entry links to the JSON definition, describes the biological motivation, and lists the goal, inputs, expected output shape, ground truth, tolerances, and the key behaviour under test.

---

## Quality Control (`atlasxomics/qc`)

### QC retention summary

- **Eval**: [filtered_cells_percentage](qc/percentage_filter.json)
- **What & Why**: Measures how many cells survive QC so downstream analyses retain enough biological diversity.
- **Goal**: Report total cells, retained cells, and retained percentage after applying QC filters.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_Healthy.h5ad`
- **Expected Output**:

  ```json
  {
    "total_cells": <int>,
    "remaining_cells": <int>,
    "remaining_percentage": <float>
  }
  ```

- **Ground Truth**: 11 027 total; 90 % retained.
- **Tolerance**: ±2 000 for totals; ±10 percentage points.
- **Core Test**: Confirms the agent computes QC metrics and applies thresholds that keep roughly 90 % of cells.

### TSSE lower tail

- **Eval**: [spatial_atac_qc_tsse_p10_numeric_tolerance_v1](qc/tsse_p10_numeric_tolerance.json)
- **What & Why**: TSS enrichment quantifies promoter accessibility; a healthy dataset retains a strong low tail.
- **Goal**: Compute the 10th percentile TSSE after QC filtering.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_Healthy.h5ad`
- **Expected Output**: `{ "p10_tsse": <float> }`.
- **Ground Truth**: 2.0.
- **Tolerance**: Minimum threshold (must be ≥ 2.0).
- **Core Test**: Ensures the agent filters out cells with very poor promoter signal.

### FRiP lower tail

- **Eval**: [spatial_atac_qc_frip_numeric_tolerance_v1](qc/frip.json)
- **What & Why**: Fraction of reads in peaks captures signal-to-noise; low FRiP cells should be dropped.
- **Goal**: Report the 10th percentile FRiP after QC.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_Healthy.h5ad`
- **Expected Output**: `{ "p10_frip": <float> }`.
- **Ground Truth**: 0.30.
- **Tolerance**: Minimum threshold (≥ 0.30).
- **Core Test**: Validates that QC preserves a strong peak signal across the cell population.

### Mitochondrial fraction

- **Eval**: [spatial_atac_qc_mito_fraction_numeric_tolerance_v1](qc/mito_numeric_tolerance.json)
- **What & Why**: High mitochondrial fractions flag damaged cells; mean levels should remain low.
- **Goal**: Compute the mean mitochondrial fraction of retained cells.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_Healthy.h5ad`
- **Expected Output**: `{ "mean_mito_fraction": <float> }`.
- **Ground Truth**: 0.20.
- **Tolerance**: Maximum threshold (≤ 0.20).
- **Core Test**: Confirms the agent applies sensible mito filters and summarises the result.

### Fragment depth per cell

- **Eval**: [spatial_atac_qc_n_fragments_numeric_tolerance_v1](qc/peaks_count.json)
- **What & Why**: Adequate fragment depth is necessary for meaningful peak calling and clustering.
- **Goal**: Return total fragments and mean `n_fragments` after QC.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_Healthy.h5ad`
- **Expected Output**:

  ```json
  {
    "total_fragments": <int>,
    "mean_n_fragments": <float>
  }
  ```

- **Ground Truth**: Mean = 10 000.
- **Tolerance**: Minimum threshold (≥ 10 000 for mean).
- **Core Test**: Checks that low-depth cells are removed.

### Fragment category balance

- **Eval**: [spatial_atac_qc_fragment_categories_distribution_v1](qc/fragment_categories_distribution.json)
- **What & Why**: Fragment length spectra (nucleosome-free, mono-, di-, tri-nucleosome) diagnose Tn5 digestion quality.
- **Goal**: Report fragment percentages per category plus the total fragment count.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_Healthy.h5ad`
- **Expected Output**: JSON with `total_cells` (total fragments) and `cell_type_distribution` of four categories.
- **Ground Truth**: 191 280 035 total fragments; percentages 60 %, 25 %, 10 %, 5 %.
- **Tolerance**: ±10 000 000 fragments; ±10 percentage points per category.
- **Core Test**: Ensures the agent reads fragment matrices and reproduces expected digestion patterns.

---

## Clustering (`atlasxomics/clustering`)

### Coarse cluster count

- **Eval**: [atlasxomics_coarse_clustering_threshold_v1](clustering/coarse_threshold.json)
- **What & Why**: Tests whether the agent can generate a low-resolution clustering that captures broad tissue compartments.
- **Goal**: Run dimensionality reduction + Leiden/Louvain at low resolution and report number of clusters discovered.
- **Input**: `latch://38438.account/AtlasxOmics/Human_DRG_Price_plotsLite/combined_sm_ge.h5ad`
- **Expected Output**: `{ "num_clusters": <int> }`.
- **Ground Truth**: 6 clusters.
- **Tolerance**: ±1 cluster.
- **Core Test**: Validates coarse clustering behaviour and reporting discipline.

### Cluster × sample independence

- **Eval**: [atlasxomics_cluster_condition_contingency_v1](clustering/batch.json)
- **What & Why**: Evaluates whether clusters are dominated by specific samples; low AMI implies clustering reflects biology over batch.
- **Goal**: Construct a cluster × sample contingency table and keep adjusted mutual information ≤ 0.2.
- **Input**: `latch://38438.account/AtlasxOmics/Human_DRG_Price_plotsLite/combined_sm_ge.h5ad`
- **Expected Output**:

  ```json
  {
    "clusters": [...],
    "samples": [...],
    "contingency": [[...], [...]]
  }
  ```

- **Ground Truth**: AMI computed from the table must satisfy ami_max = 0.2.
- **Tolerance**: Maximum threshold (≤ 0.2 AMI).
- **Core Test**: Confirms clustering does not simply recapitulate sample batching.

---

## Cell Type Annotation & Validation (`atlasxomics/cell_typing`)

### Brain major-class mapping (markers only)

- **Eval**: [mouse_brain_major_celltypes_marker_only_v2](../cell_typing/mouse_brain_major_celltypes_marker_only_v2.json)
- **What & Why**: Ensures the agent can recover the six major mouse brain cell classes using only expression-based marker heuristics.
- **Goal**: Map annotations onto the six-class vocabulary and report unique labels.
- **Input**: `latch://38438.account/Mouse_Aging_subsetted.h5ad`
- **Expected Output**: `{ "cell_types_predicted": [...] }`.
- **Ground Truth**: Immune, Neuron, Astrocyte, Endothelial, Oligodendrocyte, Microglia.
- **Tolerance**: Jaccard ≥ 0.90.
- **Core Test**: Verifies that marker scoring reproduces expected major classes without extraneous types.

### Brain vocabulary alignment (ranked genes)

- **Eval**: [brain_celltype_annotation_major_map](cell_typing/brain_celltype_annotation_major_map.json)
- **What & Why**: Confirms the agent can interpret ranked differential tables and project them onto a curated ontology.
- **Goal**: Produce mapped labels limited to the provided vocabulary.
- **Input**: Ranked genes CSV and downsampled AnnData specified in JSON.
- **Expected Output**: `{ "cell_types_predicted": [...] }` with mapped labels.
- **Ground Truth**: Same six major classes.
- **Tolerance**: Jaccard ≥ 0.90.
- **Core Test**: Exercises ontology alignment from DE results.

### Per-condition cell-type consistency

- **Eval**: [celltype_proportion_consistency_tbi_24hr_v1](cell_typing/cell_type_consistency_per_condition.json)
- **What & Why**: Tests whether annotated proportions are stable across TBI_24hr replicates, avoiding replicate-specific drift.
- **Goal**: Subset to TBI_24hr, compute per-sample cell-type proportions, and return structured JSON.
- **Input**: `latch://38438.account/AtlasxOmics/Pieper_154_brain_ArchRFull_10Core/annotated/TBI_24hr.annotated.h5ad`
- **Expected Output**: As declared in the task (`cell_types` list and sample proportion objects).
- **Ground Truth**: Chi-square q-value ≥ 0.05 and Cramér’s V ≤ 0.15.
- **Tolerance**: q-value must be ≥ 0.05 and Cramér’s V must be ≤ 0.15.
- **Core Test**: Ensures cell-type assignments are consistent across replicates of the same condition.

### Kidney 20-population composition

- **Eval**: [kidney_major_20populations_distribution](../kidney_major_20populations_distribution.json)
- **What & Why**: Checks the agent’s ability to quantify a detailed kidney atlas containing 20 populations.
- **Goal**: Report total cells and percentage for each population in the provided vocabulary.
- **Input**: `latch://38438.account/Scratch/xenium/Output/xenium_combined_12samples.h5ad`
- **Expected Output**: `{ "total_cells": <int>, "cell_type_distribution": {"TAL": <float>, ...} }`.
- **Ground Truth**: Percentages supplied in JSON (e.g., TAL 14.8091 %).
- **Tolerance**: ±100 cells; ±3 percentage points per population.
- **Core Test**: Validates high-resolution compositional analysis.

### Hepatic stellate cluster retrieval

- **Eval**: [human_liver_hepatic_stellate_clusters_v1](cell_typing/hsc_typing_1.json)
- **What & Why**: Ensures the agent can pinpoint hepatic stellate clusters within a human liver ATAC dataset.
- **Goal**: Return cluster labels that correspond to hepatic stellate cells.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_sm_ge.h5ad`
- **Expected Output**: `{ "cell_types_predicted": ["5", ...] }`.
- **Ground Truth**: Cluster “5”.
- **Tolerance**: Jaccard threshold 1.0.
- **Core Test**: Confirms lineage-specific cluster identification.

### Top-20 marker recovery panels

- **Evals**: `mouse_brain_markers_*.json`

  - [Endothelial](cell_typing/mouse_brain_markers_endothelial.json)
  - [GABAergic interneuron](cell_typing/mouse_brain_markers_gabaergic.json)
  - [Glutamatergic neuron](cell_typing/mouse_brain_markers_glutamatergic.json)
  - [Microglia](cell_typing/mouse_brain_markers_microglia.json)
  - [Neural progenitor](cell_typing/mouse_brain_markers_neural_progenitor.json)
  - [Oligodendrocyte](cell_typing/mouse_brain_markers_oligodendrocyte.json)

- **What & Why**: Each test validates the agent’s ability to recover hallmark markers for the specified lineage via 1-vs-all differential analysis.
- **Goal**: Return an ordered list `top_marker_genes` (length 20) ranked by log fold-change + Wilcoxon score.
- **Input**: `latch://38438.account/AtlasxOmics/Pieper_154_brain_ArchRFull_10Core/annotated/D01900_NG05698_TBI_24hr.annotated.h5ad`
- **Expected Output**: `{ "top_marker_genes": ["Gene1", …, "Gene20"] }`.
- **Ground Truth**: Canonical marker sets defined in each JSON.
- **Tolerance**: Precision@20 must be ≥ 0.60 and Recall@20 must be ≥ 0.50.
- **Core Test**: Ensures DE ranking surfaces known lineage markers.

### Canonical marker separation panels

- **Evals**: `*_canonical_expression_separation_v1.json`

  - [Endothelial](cell_typing/mouse_brain_endothelial_canonical_expression_separation_v1.json)
  - [GABAergic interneuron](cell_typing/mouse_brain_gabaergic_canonical_expression_separation_v1.json)
  - [Glutamatergic neuron](cell_typing/mouse_brain_glutamatergic_canonical_expression_separation_v1.json)
  - [Microglia](cell_typing/mouse_brain_microglia_canonical_expression_separation_v1.json)
  - [Neural progenitor](cell_typing/mouse_brain_neural_progenitor_canonical_expression_separation_v1.json)
  - [Oligodendrocyte](cell_typing/mouse_brain_oligodendrocyte_canonical_expression_separation_v1.json)

- **What & Why**: Quantifies the AUROC of canonical markers for each target lineage to confirm clean separation from other clusters.
- **Goal**: Return per-gene AUROC values and the mean AUROC across markers.
- **Input**: `latch://38438.account/AtlasxOmics/Pieper_154_brain_ArchRFull_10Core/annotated/D01900_NG05698_TBI_24hr.annotated.h5ad`
- **Expected Output**: `{ "per_gene_stats": [{"gene": ..., "auroc": ...}, ...], "mean_auroc": <float> }`.
- **Ground Truth**: Mean AUROC must be ≥ 0.85 and ≥ 70 % of genes must have AUROC ≥ 0.80.
- **Tolerance**: Mean AUROC must be ≥ 0.85 and at least 70 % of markers must reach AUROC ≥ 0.80.
- **Core Test**: Demonstrates strong marker separation for lineage-specific clusters.

---

## Differential Accessibility (`atlasxomics/de`)
