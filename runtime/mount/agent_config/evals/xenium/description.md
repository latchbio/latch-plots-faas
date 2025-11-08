# Xenium Spatial Single-Cell Transcriptomics Workflow + Evals

This document summarises the spatial single-cell RNA-seq evaluations that live inside:

- `xenium/qc`
- `xenium/preprocessing`
- `xenium/cell_typing`

For each evaluation you will find the biological motivation, the analytic goal, input datasets, expected output format, reference values, pass/fail tolerances, and the behaviour being tested.

---

## Quality Control (`xenium/qc`)

### Genes per cell

- **Eval**: [xenium_genes_per_cell](qc/xenium_genes_per_cell.json)
- **What & Why**: The number of detected genes per cell is a fundamental quality metric in spatial transcriptomics. Xenium datasets should show adequate gene detection across on-tissue cells, with typical values around 40-50 genes per cell for targeted panels. Monitoring these statistics ensures the agent correctly identifies on-tissue cells and applies appropriate QC thresholds.
- **Goal**: Calculate the number of detected genes (genes with at least 1 count) for each on-tissue bead and report summary statistics: total on-tissue cells, mean, median, and standard deviation of genes per cell.
- **Input**: `latch://38438.account/curio_ovary/adata_after_background_removal.h5ad`
- **Expected Output**:

  ```json
  {
    "on_tissue_cells": 1408604,
    "mean_genes_per_cell": 44.6,
    "median_genes_per_cell": 44.0,
    "std_genes_per_cell": 15.0
  }
  ```

- **Ground Truth**: On-tissue cells = 1,408,604; mean = 44.6; median = 44.0; std = 15.0.
- **Tolerance**: On-tissue cells within ±500; mean/median/std within ±5.0.
- **Core Test**: Confirms the agent correctly identifies on-tissue cells and computes gene detection statistics that match expected Xenium panel performance.

### Total UMI counts per cell

- **Eval**: [xenium_total_umi_counts](qc/xenium_total_umi_counts.json)
- **What & Why**: Total UMI counts reflect sequencing depth and cell quality. Xenium panels typically yield lower total counts per cell compared to full-transcriptome assays, but adequate depth is essential for reliable cell-type identification. Tracking UMI statistics helps verify that the dataset meets quality standards for downstream analysis.
- **Goal**: Calculate the total UMI counts (sum of all gene counts) for each cell and report summary statistics: total on-tissue cells, mean, median, and standard deviation of total counts per cell.
- **Input**: `latch://38438.account/Scratch/xenium/Output/xenium_combined_12samples.h5ad`
- **Expected Output**:

  ```json
  {
    "on_tissue_cells": 1405338,
    "mean_total_counts": 139.5,
    "median_total_counts": 118.0,
    "std_total_counts": 92.9
  }
  ```

- **Ground Truth**: On-tissue cells = 1,405,338; mean = 139.5; median = 118.0; std = 92.9.
- **Tolerance**: On-tissue cells within ±500; mean/median within ±15.0; std within ±30.0.
- **Core Test**: Validates that the agent correctly computes UMI count statistics consistent with Xenium panel sequencing depth.

---

## Preprocessing (`xenium/preprocessing`)

### Normalization and log transformation

- **Eval**: [xenium_normalization](preprocessing/xenium_normalization.json)
- **What & Why**: Normalization is critical for removing technical variation in sequencing depth across cells. Total count normalization followed by log transformation is standard for single-cell RNA-seq and ensures that downstream analyses (PCA, clustering) are not dominated by highly expressed genes or sequencing depth differences.
- **Goal**: Apply total count normalization (target_sum=1e4) and log1p transformation to the QC-filtered dataset. Report the mean sum per cell after transformations, and the min/max values in the expression matrix.
- **Input**: `latch://38438.account/curio_ovary/adata_after_qc_filtering.h5ad`
- **Expected Output**:

  ```json
  {
    "mean_sum_per_cell": 238.95,
    "max_value": 8.33,
    "min_value": 0.0
  }
  ```

- **Ground Truth**: Mean sum per cell = 238.95; max = 8.33; min = 0.0.
- **Tolerance**: Mean sum within ±50.0; max within ±1.5; min = 0.0.
- **Core Test**: Confirms the agent applies standard normalization and log transformation correctly, producing expected value ranges.

### Principal component analysis

- **Eval**: [xenium_pca](preprocessing/xenium_pca.json)
- **What & Why**: PCA reduces dimensionality while preserving major sources of variation in gene expression. For Xenium datasets with targeted gene panels (~300 genes), PCA captures cell-type-specific expression patterns and is essential for downstream batch correction and clustering. The variance explained by the top PCs indicates how well the reduced space captures biological signal.
- **Goal**: Compute PCA on the filtered data with n_comps=50. Report the number of observations, genes, principal components computed, and the sum of variance ratios for the first 10 PCs.
- **Input**: `latch://38438.account/Scratch/xenium/evals/xenium_kidney_filtered_norm.h5ad`
- **Expected Output**:

  ```json
  {
    "n_obs": 1405315,
    "n_vars": 300,
    "n_pcs": 50,
    "variance_explained_top10": 0.3155
  }
  ```

- **Ground Truth**: n_obs = 1,405,315; n_vars = 300; n_pcs = 50; variance_explained_top10 = 0.3155.
- **Tolerance**: n_obs and n_vars must match exactly; n_pcs within ±40; variance_explained_top10 within ±0.025.
- **Core Test**: Verifies the agent computes PCA correctly and captures expected variance in the top components.

### Batch correction with Harmony

- **Eval**: [xenium_batch_correction_same_batch_fraction](preprocessing/xenium_batch_correction_same_batch_fraction.json)
- **What & Why**: Multi-sample Xenium experiments often exhibit batch effects that can confound biological signal. Harmony integrates cells across batches by correcting for technical variation while preserving biological structure. The same-batch fraction in the k-nearest neighbor graph measures integration quality: lower values indicate better mixing across batches.
- **Goal**: Given PCA embeddings and batch labels, run Harmony to produce corrected embeddings. Build a 30-NN graph on the corrected space and compute the mean fraction of neighbors from the same batch.
- **Input**: `latch://38438.account/Scratch/xenium/Output/xenium_mouse_kidney_no_celltype_postprocess_100k.h5ad`
- **Expected Output**:

  ```json
  {
    "n_batches": 12,
    "n_pcs": 50,
    "same_batch_fraction_30nn_mean": 0.1384
  }
  ```

- **Ground Truth**: n_batches = 12; n_pcs = 50; same_batch_fraction_30nn_mean = 0.1384.
- **Tolerance**: n_batches must match exactly; n_pcs within ±40; same_batch_fraction within ±0.1.
- **Core Test**: Demonstrates that Harmony successfully integrates batches, reducing same-batch clustering in the k-NN graph.

### Leiden clustering on corrected PCs

- **Eval**: [xenium_leiden_on_corrected_pcs](preprocessing/xenium_leiden_on_corrected_pcs.json)
- **What & Why**: Clustering identifies cell populations based on expression similarity. Leiden clustering on batch-corrected PCA space should reveal biologically meaningful groups without batch-driven artifacts. The cluster size distribution (entropy) and largest cluster fraction indicate whether clustering resolution is appropriate.
- **Goal**: Using the corrected 50-PC embedding, build a 30-NN graph and run Leiden clustering with resolution=1.2. Report the number of clusters, largest cluster fraction, and cluster size entropy.
- **Input**: `latch://38438.account/Scratch/xenium/Output/xenium_mouse_kidney_no_celltype_postprocess_100k.h5ad`
- **Expected Output**:

  ```json
  {
    "resolution_used": 1.2,
    "n_clusters": 15,
    "largest_cluster_frac": 0.1626,
    "cluster_size_entropy": 2.5360
  }
  ```

- **Ground Truth**: resolution = 1.2; n_clusters = 15; largest_cluster_frac = 0.1626; cluster_size_entropy = 2.5360.
- **Tolerance**: resolution must match exactly; n_clusters within ±3; largest_cluster_frac within ±0.02; entropy within ±0.2.
- **Core Test**: Confirms the agent applies Leiden clustering correctly and produces a balanced cluster distribution at the specified resolution.

---

## Cell Type Annotation & Validation (`xenium/cell_typing`)

### Major 20-population distribution

- **Eval**: [kidney_major_20populations_distribution](cell_typing/kidney_major_20populations_distribution.json)
- **What & Why**: The mouse kidney contains a rich diversity of cell types across nephron segments, vasculature, and stroma. Accurately quantifying all 20 major populations (Podocytes, Glomerular EC, Proximal Tubule segments, Loop of Henle, Collecting Duct, etc.) validates that the agent can handle complex tissue composition and correctly assign cell-type labels.
- **Goal**: Assign each cell to one of 20 populations using the provided vocabulary and report the overall distribution as percentages of total cells per type.
- **Input**: `latch://38438.account/Scratch/xenium/Output/xenium_mouse_kidney_no_celltype_postprocess_200k.h5ad`
- **Expected Output** (excerpt):

  ```json
  {
    "total_cells": 1374915,
    "cell_type_distribution": {
      "TAL": 14.8091,
      "Fib": 11.6609,
      "Inj-PT": 10.6259,
      "PTS2": 10.6108,
      "non-glom EC": 9.4014,
      "PTS1": 8.8904,
      "PTS3": 5.7096,
      "Immune": 4.5385,
      "Per-SMC": 3.9770,
      "DCT": 3.5474,
      "CNT": 2.7517,
      "PC": 2.5505,
      "Glom-EC": 1.9837,
      "ICB": 1.8797,
      "DTL": 1.7714,
      "FR-PT": 1.5767,
      "Pod": 1.2409,
      "ICA": 1.0974,
      "Uro": 0.8110,
      "PEC": 0.5661
    }
  }
  ```

- **Ground Truth**: Values shown above (total = 1,374,915 cells).
- **Tolerance**: Total cells within ±100; each percentage within ±3 percentage points.
- **Core Test**: Demonstrates the agent can reproduce fine-grained kidney composition across all major populations.

### Major 20-population label set

- **Eval**: [kidney_major_20populations_labelset](cell_typing/kidney_major_20populations_labelset.json)
- **What & Why**: Beyond distribution accuracy, the agent must identify all 20 cell types present in the dataset. This label set validation ensures complete coverage of the kidney cell-type vocabulary without missing or incorrectly mapping populations.
- **Goal**: Assign cells to the 20-population vocabulary and report the set of unique cell types present in predictions.
- **Input**: `latch://38438.account/Scratch/xenium/Output/xenium_mouse_kidney_no_celltype_postprocess_200k.h5ad`
- **Expected Output**:

  ```json
  {
    "cell_types_predicted": [
      "Pod", "Glom-EC", "EC", "PTS1", "PTS2", "PTS3",
      "Inj_PT", "FR_PT", "DTL", "TAL", "DCT", "CNT",
      "PC", "ICA", "ICB", "Uro", "PEC", "Fib",
      "Per-SMC", "Immune"
    ]
  }
  ```

- **Ground Truth**: All 20 categories listed above.
- **Tolerance**: Jaccard similarity must equal 1.0 (exact match).
- **Core Test**: Validates that the agent identifies all 20 kidney cell types without omissions or incorrect mappings.

### Proximal tubule subclustering

- **Eval**: [classify_pt_distribution_subcluster](cell_typing/classify_pt_distribution_subcluster.json)
- **What & Why**: Proximal tubule (PT) cells are heterogeneous, with distinct segments (PTS1, PTS2, PTS3) and injury states (Inj-PT, FR-PT). Subclustering PT cells and assigning segment-specific labels using marker genes (e.g., Slc5a12 for PTS1, Cyp2e1 for PTS2) tests the agent's ability to perform fine-grained subtyping within a major cell class.
- **Goal**: Cluster PT cells into 5 subclusters and assign labels using marker-guided annotation. Report the distribution of PT subtypes as percentages of total PT cells.
- **Input**: `latch://38438.account/Scratch/xenium/evals/xenium_mouse_kidney_subset_PT_day2_no_celltype.h5ad`
- **Expected Output**:

  ```json
  {
    "total_cells": 80028,
    "cell_type_distribution": {
      "Inj_PT": 48.55,
      "PTS1": 42.06,
      "PTS2": 5.02,
      "PTS3": 0.90,
      "FR_PT": 3.47
    }
  }
  ```

- **Ground Truth**: Values shown above.
- **Tolerance**: Each percentage within ±3 percentage points.
- **Core Test**: Confirms the agent can subcluster and annotate PT segments using marker gene expression.

### Proximal tubule subclustering (advanced)

- **Eval**: [classify_pt_distribution_advanced](cell_typing/classify_pt_distribution_advanced.json)
- **What & Why**: Similar to the basic PT subclustering, but requires the agent to subset cells based on metadata (cell_type containing 'PT' and time='Day2') before subclustering. This tests the agent's ability to work with pre-annotated data and perform conditional subtyping.
- **Goal**: Subset to PT cells at Day2, then cluster and annotate PT subtypes. Report the distribution as percentages.
- **Input**: `latch://38438.account/Scratch/xenium/evals/evals_output/kidney_agent_annotated_20_cells.h5ad`
- **Expected Output**: Same format as above with values: Inj_PT = 48.55%, PTS1 = 42.06%, PTS2 = 5.02%, PTS3 = 0.90%, FR_PT = 3.47%.
- **Ground Truth**: Same as above.
- **Tolerance**: Each percentage within ±3 percentage points.
- **Core Test**: Validates conditional subclustering on pre-annotated datasets.

### Collecting duct subclustering

- **Eval**: [collecting_duct_subcluster_distribution](cell_typing/collecting_duct_subcluster_distribution.json)
- **What & Why**: The collecting duct contains principal cells (PC) and two types of intercalated cells (ICA, ICB) with distinct functions. Subclustering and annotating these subtypes using canonical markers (Aqp2 for PC, Atp6v1b1/Slc4a1 for ICA, Slc26a4 for ICB) tests fine-grained annotation within a spatially restricted region.
- **Goal**: Subset to collecting duct cells, cluster into 3 subclusters, and assign labels using marker genes. Report the distribution of PC, ICA, and ICB as percentages.
- **Input**: `latch://38438.account/Scratch/xenium/evals/xenium_mouse_kidney_no_celltype_postprocess.h5ad`
- **Expected Output**:

  ```json
  {
    "total_cells": 36441,
    "cell_type_distribution": {
      "PC": 30.22,
      "ICA": 14.22,
      "ICB": 55.56
    }
  }
  ```

- **Ground Truth**: Values shown above.
- **Tolerance**: Each percentage within ±3 percentage points.
- **Core Test**: Demonstrates subclustering and marker-guided annotation within the collecting duct niche.

### Interstitial cell classification

- **Eval**: [classify_interstitial_cells_distribution](cell_typing/classify_interstitial_cells_distribution.json)
- **What & Why**: Interstitial cells in the kidney include fibroblasts (Fib) and pericytes/smooth muscle cells (Per-SMC), which play distinct roles in tissue structure and fibrosis. Correctly distinguishing these populations is important for understanding kidney pathology and spatial organization.
- **Goal**: Classify interstitial cells into Fib vs Per-SMC and report the distribution as percentages.
- **Input**: `latch://38438.account/Scratch/xenium/evals/xenium_kidney_interstitial_cells.h5ad`
- **Expected Output**:

  ```json
  {
    "total_cells": 215009,
    "cell_type_distribution": {
      "Fib": 74.57,
      "Per-SMC": 25.43
    }
  }
  ```

- **Ground Truth**: Fib = 74.57%, Per-SMC = 25.43%.
- **Tolerance**: Each percentage within ±3 percentage points.
- **Core Test**: Validates the agent can distinguish fibroblast and pericyte populations in the interstitium.

### Interstitial cell label set

- **Eval**: [classify_interstitial_cells_distribution_labelset](cell_typing/classify_interstitial_cells_distribution_labelset.json)
- **What & Why**: Ensures the agent identifies both Fib and Per-SMC populations without missing either category.
- **Goal**: Classify interstitial cells and report the set of unique cell types present.
- **Input**: `latch://38438.account/Scratch/xenium/evals/xenium_kidney_interstitial_cells.h5ad`
- **Expected Output**:

  ```json
  {
    "cell_types_predicted": ["Fib", "Per-SMC"]
  }
  ```

- **Ground Truth**: Both "Fib" and "Per-SMC".
- **Tolerance**: Jaccard similarity must equal 1.0.
- **Core Test**: Confirms complete coverage of interstitial cell types.

### PEC marker gene discovery

- **Eval**: [kidney_pec_marker_gene_top20](cell_typing/kidney_pec_marker_gene_top20.json)
- **What & Why**: Parietal epithelial cells (PEC) are specialized cells in the glomerulus. Identifying their top marker genes through differential expression analysis validates that the agent can discover cell-type-specific markers and rank them appropriately.
- **Goal**: For PEC cells, perform 1-vs-all Wilcoxon differential expression, rank by log fold change and Wilcoxon score, and report the top 20 marker genes.
- **Input**: `latch://38438.account/Scratch/xenium/evals/evals_output/kidney_agent_annotated_20_cells.h5ad`
- **Expected Output**:

  ```json
  {
    "top_marker_genes": ["Wnt16", "Akap12", "Kcnmb2", ...]
  }
  ```

- **Ground Truth**: Canonical markers include Wnt16, Akap12, Kcnmb2, C3, Kcnma1, Adgrv1, Amph, Proser2, Tgfb2, Sorcs3, Adamts16, Nphs2, Col18a1, Upk1b, Vcam1, Trpm6, Synpo, Myom2, Efnb1, Jag1.
- **Tolerance**: Precision@20 ≥ 0.60 and Recall@20 ≥ 0.50.
- **Core Test**: Confirms the agent discovers canonical PEC markers through differential expression analysis.

### Urothelial marker gene discovery

- **Eval**: [kidney_uro_marker_gene_top20](cell_typing/kidney_uro_marker_gene_top20.json)
- **What & Why**: Urothelial cells (Uro) line the renal pelvis and ureter. Identifying their markers (e.g., Krt15, Upk1b, Msln) tests marker discovery for a rare but distinct kidney population.
- **Goal**: For Uro cells, perform 1-vs-all Wilcoxon differential expression and report the top 20 marker genes.
- **Input**: `latch://38438.account/Scratch/xenium/evals/evals_output/kidney_agent_annotated_20_cells.h5ad`
- **Expected Output**:

  ```json
  {
    "top_marker_genes": ["Krt15", "Upk1b", "Msln", ...]
  }
  ```

- **Ground Truth**: Canonical markers include Krt15, Upk1b, Msln, Krt5, Krt14, Krt19, Abcc3, Ano1, Klf5, Efnb2, Akr1b3, Runx1, Gprc5a, Spp1, Muc20, Kcnma1, Cd24a, Tpm1, Proser2, Tnfrsf11a.
- **Tolerance**: Precision@20 ≥ 0.60 and Recall@20 ≥ 0.50.
- **Core Test**: Validates marker discovery for urothelial cells.

### DTL/TAL marker discovery

- **Eval**: [dtl_tal_markers_wilcoxon_labelset](cell_typing/dtl_tal_markers_wilcoxon_labelset.json)
- **What & Why**: The Loop of Henle contains descending (DTL) and ascending (TAL) limbs with distinct functions and marker genes. Discovering markers for each segment (e.g., UMOD/SLC12A1 for TAL, FST/AQP1/BST1 for DTL) tests the agent's ability to perform differential expression within a related cell lineage.
- **Goal**: Identify TAL and DTL cells, perform 1-vs-rest Wilcoxon tests for each class, and return top markers ranked by log fold change and adjusted p-value.
- **Input**: `latch://38438.account/Scratch/xenium/Output/xenium_mouse_kidney_no_celltype_postprocess_200k.h5ad`
- **Expected Output**:

  ```json
  {
    "markers_predicted": {
      "TAL": ["UMOD", "SLC12A1", ...],
      "DTL": ["FST", "AQP1", "BST1", ...]
    }
  }
  ```

- **Ground Truth**: TAL markers include UMOD, SLC12A1; DTL markers include FST, AQP1, BST1.
- **Tolerance**: Jaccard similarity ≥ 0.5 for each label's top-15 markers.
- **Core Test**: Confirms the agent discovers segment-specific markers for Loop of Henle cell types.

### Glomerular vs non-glomerular EC markers

- **Eval**: [glom_ec_vs_nonglom_ec_markers_wilcoxon_labelset](cell_typing/glom_ec_vs_nonglom_ec_markers_wilcoxon_labelset.json)
- **What & Why**: Endothelial cells in the kidney are heterogeneous: glomerular EC (Glom-EC) have specialized fenestrated morphology, while non-glomerular EC (EC) serve different vascular functions. Performing 1-vs-1 differential expression between these subtypes tests the agent's ability to identify subtle distinctions within a cell lineage.
- **Goal**: Split endothelial cells into Glom-EC vs EC, perform 1-vs-1 Wilcoxon tests, and return top markers for each class.
- **Input**: `latch://38438.account/Scratch/xenium/Output/xenium_mouse_kidney_no_celltype_postprocess_200k.h5ad`
- **Expected Output**:

  ```json
  {
    "markers_predicted": {
      "Glom-EC": ["EHD3", "PLVAP", ...],
      "EC": ["VWF", "ICAM2", ...]
    }
  }
  ```

- **Ground Truth**: Glom-EC markers include EHD3, PLVAP; EC markers include VWF, ICAM2.
- **Tolerance**: Jaccard similarity ≥ 0.5 for each label's top-15 markers.
- **Core Test**: Validates marker discovery for endothelial cell subtypes using 1-vs-1 comparison.

### Podocyte marker expression separation

- **Eval**: [kidney_podocyte_nphs2_expression_separation_v1](cell_typing/kidney_podocyte_nphs2_expression_separation_v1.json)
- **What & Why**: Nphs2 (nephrin) is a canonical podocyte marker essential for glomerular filtration. Computing AUROC for Nphs2 expression in podocytes vs all other cell types measures how cleanly this marker distinguishes podocytes, validating annotation quality.
- **Goal**: Subset to Hour4 timepoint, compute AUROC for Nphs2 expression in Podocyte vs all other cell types.
- **Input**: `latch://38438.account/Scratch/xenium/evals/evals_output/kidney_agent_annotated_20_cells.h5ad`
- **Expected Output**:

  ```json
  {
    "per_gene_stats": [
      {"gene": "Nphs2", "auroc": 0.95}
    ],
    "mean_auroc": 0.95
  }
  ```

- **Ground Truth**: AUROC(Nphs2: Podocyte vs rest) ≥ 0.90.
- **Tolerance**: Mean AUROC must be ≥ 0.90; with one gene, fraction_high = 1.0 enforces that Nphs2 clears the 0.90 cutoff.
- **Core Test**: Ensures canonical podocyte markers show strong expression separation in annotated podocytes.

### Fibroblast marker expression separation

- **Eval**: [kidney_fibroblasts_canonical_expression_separation](cell_typing/kidney_fibroblasts_canonical_expression_separation.json)
- **What & Why**: Col1a1 (collagen type I alpha 1) is a canonical fibroblast marker. Computing AUROC for Col1a1 in fibroblasts vs all other cell types validates that fibroblast annotations align with expected marker expression patterns.
- **Goal**: Compute AUROC for Col1a1 expression in Fib vs all other cell types.
- **Input**: `latch://38438.account/Scratch/xenium/evals/evals_output/kidney_agent_annotated_20_cells.h5ad`
- **Expected Output**:

  ```json
  {
    "per_gene_stats": [
      {"gene": "Col1a1", "auroc": 0.88}
    ],
    "mean_auroc": 0.88
  }
  ```

- **Ground Truth**: Mean AUROC ≥ 0.85; ≥70% of markers must have AUROC ≥ 0.80 (single marker case).
- **Tolerance**: Mean AUROC must be ≥ 0.85 and ≥70% of markers must reach AUROC ≥ 0.80.
- **Core Test**: Confirms canonical fibroblast markers show strong expression separation.

### Pericyte marker expression separation

- **Eval**: [kidney_pericytes_canonical_expression_separation](cell_typing/kidney_pericytes_canonical_expression_separation.json)
- **What & Why**: Myh11 (myosin heavy chain 11) is a canonical marker for pericytes and smooth muscle cells. Validating expression separation for this marker ensures that Per-SMC annotations are biologically meaningful.
- **Goal**: Compute AUROC for Myh11 expression in Per-SMC vs all other cell types.
- **Input**: `latch://38438.account/Scratch/xenium/evals/evals_output/kidney_agent_annotated_20_cells.h5ad`
- **Expected Output**: Same format as above with mean AUROC ≥ 0.85.
- **Ground Truth**: Mean AUROC ≥ 0.85; ≥70% of markers must have AUROC ≥ 0.80.
- **Tolerance**: Mean AUROC must be ≥ 0.85 and ≥70% of markers must reach AUROC ≥ 0.80.
- **Core Test**: Validates pericyte marker expression separation.

### IC/PC spatial adjacency

- **Eval**: [ic_pc_spatial_adjacency](cell_typing/ic_pc_spatial_adjacency.json)
- **What & Why**: Intercalated cells (IC) and principal cells (PC) are intermingled in the collecting duct, with IC cells typically located within 15-55 µm of PC cells. Quantifying spatial adjacency validates that annotations respect known tissue architecture and tests the agent's ability to perform spatial analysis.
- **Goal**: For each IC cell, find its nearest PC and compute IC→PC distances. Report median, 90th percentile distances, and the percentage of IC cells within 15 µm and 55 µm of a PC.
- **Input**: `latch://38438.account/Scratch/xenium/evals/evals_output/kidney_agent_annotated_20_cells.h5ad`
- **Expected Output**:

  ```json
  {
    "median_ic_to_pc_um": 12.5,
    "p90_ic_to_pc_um": 45.0,
    "pct_ic_within_15um": 75.0,
    "pct_ic_mixed_within_55um": 90.0,
    "adjacency_pass": true
  }
  ```

- **Ground Truth**: Median ≤ 25 µm, 90th percentile ≤ 80 µm, ≥60% of IC within 15 µm, ≥60% of IC within 55 µm.
- **Tolerance**: Pass if all thresholds are met.
- **Core Test**: Confirms IC and PC cells show expected spatial proximity, validating annotation and spatial analysis capabilities.

### FR-PT composition over time

- **Eval**: [fr_pt_composition_over_time](cell_typing/fr_pt_composition_over_time.json)
- **What & Why**: Failed repair proximal tubule (FR-PT) cells are a pathological state that increases after kidney injury. Tracking FR-PT fraction over time (Sham → 4h → 12h → 2d → 14d → 6w) tests the agent's ability to perform temporal analysis and detect biologically meaningful trends in cell-type composition.
- **Goal**: Compute FR-PT fraction per time group and determine if the fraction increases by 14d compared to baseline (Sham). Report per-group counts and percentages.
- **Input**: `latch://38438.account/Scratch/xenium/evals/evals_output/kidney_agent_annotated_20_cells.h5ad`
- **Expected Output**:

  ```json
  {
    "per_group": [
      {"group": "Sham", "count_target": 100, "count_all": 5000, "percent": 2.0},
      {"group": "4h", "count_target": 150, "count_all": 5000, "percent": 3.0},
      ...
    ],
    "increases_by_compare_group": true,
    "delta_vs_baseline_pct": 5.0
  }
  ```

- **Ground Truth**: FR-PT fraction at 14d is significantly higher than baseline (Fisher exact test, one-sided, α=0.05).
- **Tolerance**: Pass if statistical test indicates significant increase.
- **Core Test**: Validates temporal analysis and detection of injury-related cell-type changes.

