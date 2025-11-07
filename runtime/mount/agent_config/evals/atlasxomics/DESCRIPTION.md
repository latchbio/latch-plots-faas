# AtlasXomics Spatial ATAC-seq Workflow + Evals

This document summarises the spatial single-cell ATAC-seq evaluations that live inside:

- `atlasxomics/qc`
- `atlasxomics/clustering`
- `atlasxomics/cell_typing`
- `atlasxomics/de`

For each evaluation you will find the biological motivation, the analytic goal, input datasets, expected output format, reference values, pass/fail tolerances, and the behaviour being tested.

---

## Quality Control (`atlasxomics/qc`)

### QC retention summary

- **Eval**: [filtered_cells_percentage](qc/percentage_filter.json)
- **What & Why**: Spatial ATAC experiments typically retain the majority of profiled nuclei after gentle filtering. Tracking how many cells remain after QC ensures the agent applies sensible thresholds instead of over-pruning the dataset.
- **Goal**: Report the total number of cells, the number retained after QC, and the retained percentage.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_Healthy.h5ad`
- **Expected Output**:

  ```json
  {
    "total_cells": 11027,
    "remaining_cells": 9924,
    "remaining_percentage": 90.0
  }
  ```

- **Ground Truth**: Total cells = 11 027; retained cells = 9 924; retained percentage = 90.0 %.
- **Tolerance**: Total cells within ±2 000; retained percentage within ±10 percentage points.
- **Core Test**: Confirms the agent inspects QC metrics and selects thresholds that preserve ~90 % of cells.

### TSSE lower tail

- **Eval**: [spatial_atac_qc_tsse_p10_numeric_tolerance_v1](qc/tsse_p10_numeric_tolerance.json)
- **What & Why**: TSS enrichment (TSSE) reflects promoter accessibility; healthy nuclei show strong enrichment while compromised ones do not. Ensuring the 10th percentile stays high verifies that poor-quality cells are discarded.
- **Goal**: Compute the 10th percentile TSSE after QC filtering.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_Healthy.h5ad`
- **Expected Output**:

  ```json
  {
    "p10_tsse": 2.0
  }
  ```

- **Ground Truth**: p10 TSSE = 2.0.
- **Tolerance**: Must be ≥ 2.0.
- **Core Test**: Ensures the agent flags cells whose promoter signal falls below the accepted bound.

### FRiP lower tail

- **Eval**: [spatial_atac_qc_frip_numeric_tolerance_v1](qc/frip.json)
- **What & Why**: The fraction of reads in peaks (FRiP) measures signal-to-noise. Low FRiP indicates diffuse background rather than true accessibility. Monitoring the 10th percentile confirms that the majority of cells retain strong peak signal.
- **Goal**: Report the 10th percentile FRiP after QC.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_Healthy.h5ad`
- **Expected Output**:

  ```json
  {
    "p10_frip": 0.3
  }
  ```

- **Ground Truth**: p10 FRiP = 0.30.
- **Tolerance**: Must be ≥ 0.30.
- **Core Test**: Validates that noisy cells are removed while peak-rich cells remain.

### Mitochondrial fraction

- **Eval**: [spatial_atac_qc_mito_fraction_numeric_tolerance_v1](qc/mito_numeric_tolerance.json)
- **What & Why**: Damaged or leaking nuclei exhibit elevated mitochondrial DNA accessibility. Controlling the mean mitochondrial fraction keeps the dataset focused on healthy chromatin profiles.
- **Goal**: Compute the mean mitochondrial fraction of retained cells.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_Healthy.h5ad`
- **Expected Output**:

  ```json
  {
    "mean_mito_fraction": 0.2
  }
  ```

- **Ground Truth**: Mean mitochondrial fraction = 0.20.
- **Tolerance**: Must be ≤ 0.20.
- **Core Test**: Checks that mitochondrial contamination is kept within acceptable limits.

### Fragment depth per cell

- **Eval**: [spatial_atac_qc_n_fragments_numeric_tolerance_v1](qc/peaks_count.json)
- **What & Why**: Adequate fragment depth underpins reliable peak discovery and downstream clustering. Averaging fragment counts across cells confirms the dataset was sequenced deeply enough.
- **Goal**: Report the total fragment count and the mean `n_fragments` after QC.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_Healthy.h5ad`
- **Expected Output**:

  ```json
  {
    "total_fragments": 117831111,
    "mean_n_fragments": 10000.0
  }
  ```

- **Ground Truth**: Mean `n_fragments` = 10 000; total fragments = 117 831 111.
- **Tolerance**: Mean must be ≥ 10 000; total fragments within ±10 % (implicitly covered by the mean threshold).
- **Core Test**: Ensures low-depth cells are removed and sufficient fragments remain for peak calling.

### Fragment category balance

- **Eval**: [spatial_atac_qc_fragment_categories_distribution_v1](qc/fragment_categories_distribution.json)
- **What & Why**: Spatial ATAC libraries reveal characteristic fragment-length signatures (nucleosome-free vs mono/di/tri nucleosome). Deviations signal over/under digestion or technical failure.
- **Goal**: Quantify the percentage of fragments in each length category and report the total fragment count.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_Healthy.h5ad`
- **Expected Output**:

  ```json
  {
    "total_cells": 191280035,
    "cell_type_distribution": {
      "Nucleosome-free": 60.0,
      "Mono-nucleosome": 25.0,
      "Di-nucleosome": 10.0,
      "Tri-nucleosome": 5.0
    }
  }
  ```

- **Ground Truth**: Total fragments = 191 280 035; percentages 60 % nucleosome-free, 25 % mono, 10 % di, 5 % tri.
- **Tolerance**: Total within ±10 000 000; each percentage within ±10 percentage points.
- **Core Test**: Demonstrates that the agent extracts fragment-length spectra consistent with well-prepared libraries.

---

## Clustering (`atlasxomics/clustering`)

### Coarse cluster count

- **Eval**: [atlasxomics_coarse_clustering_threshold_v1](clustering/coarse_threshold.json)
- **What & Why**: High-level clustering should capture broad tissue states (e.g., dorsal vs ventral). Counting coarse clusters confirms the agent explores appropriate resolutions rather than over/under splitting the tissue.
- **Goal**: Perform dimensionality reduction, run Leiden/Louvain at low resolution, and report the number of clusters discovered.
- **Input**: `latch://38438.account/AtlasxOmics/Human_DRG_Price_plotsLite/combined_sm_ge.h5ad`
- **Expected Output**:

  ```json
  {
    "num_clusters": 6
  }
  ```

- **Ground Truth**: 6 clusters.
- **Tolerance**: ±1 cluster.
- **Core Test**: Verifies that the agent adjusts resolution parameters to produce a sensible coarse partition.

### Cluster × sample independence

- **Eval**: [atlasxomics_cluster_condition_contingency_v1](clustering/batch.json)
- **What & Why**: Biological clusters should not simply replicate sample batches. Computing adjusted mutual information (AMI) between clusters and samples checks that clustering is not driven by technical factors.
- **Goal**: Create a contingency table counting cells per cluster/sample pair and ensure AMI ≤ 0.2.
- **Input**: `latch://38438.account/AtlasxOmics/Human_DRG_Price_plotsLite/combined_sm_ge.h5ad`
- **Expected Output**:

  ```json
  {
    "clusters": ["0", "1", "2", "3", "4", "5"],
    "samples": ["Sample_A", "Sample_B", "Sample_C"],
    "contingency": [
      [812, 430, 377],
      [642, 501, 318],
      [701, 462, 299],
      [588, 544, 261],
      [455, 379, 274],
      [407, 352, 229]
    ]
  }
  ```

- **Ground Truth**: The contingency table above yields AMI = 0.18.
- **Tolerance**: AMI must be ≤ 0.20.
- **Core Test**: Confirms clusters are biologically meaningful rather than sample-specific artefacts.

---

## Cell Type Annotation & Validation (`atlasxomics/cell_typing`)

### Brain major-class mapping (markers only)

- **Eval**: [mouse_brain_major_celltypes_marker_only_v2](../cell_typing/mouse_brain_major_celltypes_marker_only_v2.json)
- **What & Why**: Marker-based annotation tests whether the agent recognises the six dominant brain lineages (Neurons, Astrocytes, Oligodendrocytes, Microglia, Endothelial cells, Immune cells) without spatial cues.
- **Goal**: Map marker-based annotations to the six-class vocabulary and return the unique labels.
- **Input**: `latch://38438.account/Mouse_Aging_subsetted.h5ad`
- **Expected Output**:

  ```json
  {
    "cell_types_predicted": [
      "Neuron",
      "Astrocyte",
      "Oligodendrocyte",
      "Microglia",
      "Endothelial",
      "Immune"
    ]
  }
  ```

- **Ground Truth**: The six classes above.
- **Tolerance**: Jaccard similarity ≥ 0.90 against the six-class set.
- **Core Test**: Demonstrates the agent can translate canonical marker panels into high-level annotations.

### Brain vocabulary alignment (ranked genes)

- **Eval**: [brain_celltype_annotation_major_map](cell_typing/brain_celltype_annotation_major_map.json)
- **What & Why**: Ranked differential lists (per sample) must be collated and mapped to a controlled ontology. This eval checks the agent’s ability to digest white-listed rank tables.
- **Goal**: Produce mapped labels restricted to the provided vocabulary and report them.
- **Input**:
  - `latch://38438.account/AtlasxOmics/Pieper_154_brain_ArchRFull_10Core/all_ranked_genes.csv`
  - `latch://38438.account/AtlasxOmics/Pieper_154_brain_ArchRFull_10Core/combined_vols_sm_ge_downsampled_100k.h5ad`
- **Expected Output**:

  ```json
  {
    "cell_types_predicted": [
      "Neuron",
      "Astrocyte",
      "Oligodendrocyte",
      "Microglia",
      "Endothelial",
      "Immune"
    ]
  }
  ```

- **Ground Truth**: Same six major classes.
- **Tolerance**: Jaccard similarity ≥ 0.90.
- **Core Test**: Validates ontology mapping when the agent is supplied with pre-ranked marker genes.

### Per-condition cell-type consistency

- **Eval**: [celltype_proportion_consistency_tbi_24hr_v1](cell_typing/cell_type_consistency_per_condition.json)
- **What & Why**: Replicates within a condition should display similar cell-type compositions. Statistical tests guard against condition-specific annotation drift.
- **Goal**: Subset to TBI_24hr, compute per-sample cell-type proportions, and return a structured summary.
- **Input**: `latch://38438.account/AtlasxOmics/Pieper_154_brain_ArchRFull_10Core/annotated/TBI_24hr.annotated.h5ad`
- **Expected Output**:

  ```json
  {
    "cell_types": ["Neuron", "Astrocyte", "Oligodendrocyte", "Microglia", "Endothelial", "Immune"],
    "samples": [
      {
        "sample": "D01876_NG05324",
        "proportions": {
          "Neuron": 0.47,
          "Astrocyte": 0.18,
          "Oligodendrocyte": 0.12,
          "Microglia": 0.09,
          "Endothelial": 0.08,
          "Immune": 0.06
        }
      },
      {
        "sample": "D01877_NG05325",
        "proportions": {
          "Neuron": 0.45,
          "Astrocyte": 0.20,
          "Oligodendrocyte": 0.13,
          "Microglia": 0.08,
          "Endothelial": 0.08,
          "Immune": 0.06
        }
      },
      {
        "sample": "D01878_NG05326",
        "proportions": {
          "Neuron": 0.46,
          "Astrocyte": 0.19,
          "Oligodendrocyte": 0.11,
          "Microglia": 0.09,
          "Endothelial": 0.09,
          "Immune": 0.06
        }
      }
    ]
  }
  ```

- **Ground Truth**: q-value = 0.32 (Benjamini–Hochberg) and Cramér’s V = 0.11.
- **Tolerance**: q-value must be ≥ 0.05 and Cramér’s V must be ≤ 0.15.
- **Core Test**: Confirms that annotated proportions are statistically consistent across replicates.

### Kidney 20-population composition

- **Eval**: [kidney_major_20populations_distribution](../kidney_major_20populations_distribution.json)
- **What & Why**: Spatial kidney datasets contain a rich variety of tubular, stromal, and immune populations. Quantifying all 20 reference populations while staying within tolerance indicates the agent can handle complex compositional analyses.
- **Goal**: Report total cell count and percentage contribution for each of the 20 populations.
- **Input**: `latch://38438.account/Scratch/xenium/Output/xenium_combined_12samples.h5ad`
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

- **Ground Truth**: Values shown above.
- **Tolerance**: Total cells within ±100; each percentage within ±3 percentage points.
- **Core Test**: Demonstrates the agent can reproduce fine-grained kidney composition.

### Hepatic stellate cluster retrieval

- **Eval**: [human_liver_hepatic_stellate_clusters_v1](cell_typing/hsc_typing_1.json)
- **What & Why**: Hepatic stellate cells regulate fibrosis in liver injury. Correctly identifying their clusters verifies lineage-specific annotation skills.
- **Goal**: Return the cluster IDs corresponding to hepatic stellate cells after gene activity scoring.
- **Input**: `latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC/combined_sm_ge.h5ad`
- **Expected Output**:

  ```json
  {
    "cell_types_predicted": ["5"]
  }
  ```

- **Ground Truth**: Cluster “5”.
- **Tolerance**: Jaccard similarity must equal 1.0.
- **Core Test**: Ensures the agent pinpoints the known hepatic stellate cluster without spurious labels.

### Top-20 marker recovery panels

- **Evals**: `mouse_brain_markers_*.json`

  - [Endothelial](cell_typing/mouse_brain_markers_endothelial.json)
  - [GABAergic interneuron](cell_typing/mouse_brain_markers_gabaergic.json)
  - [Glutamatergic neuron](cell_typing/mouse_brain_markers_glutamatergic.json)
  - [Microglia](cell_typing/mouse_brain_markers_microglia.json)
  - [Neural progenitor](cell_typing/mouse_brain_markers_neural_progenitor.json)
  - [Oligodendrocyte](cell_typing/mouse_brain_markers_oligodendrocyte.json)

- **What & Why**: Differential marker recovery demonstrates that cluster labels align with established lineage markers (e.g., PECAM1 for endothelial, GAD1 for GABAergic neurons).
- **Goal**: Run 1-vs-all Wilcoxon tests and return the top 20 markers in ranked order for each lineage.
- **Input**: `latch://38438.account/AtlasxOmics/Pieper_154_brain_ArchRFull_10Core/annotated/D01900_NG05698_TBI_24hr.annotated.h5ad`
- **Expected Output**:

  ```json
  {
    "top_marker_genes": [
      "Gene1",
      "Gene2",
      "Gene3",
      "Gene4",
      "Gene5",
      "Gene6",
      "Gene7",
      "Gene8",
      "Gene9",
      "Gene10",
      "Gene11",
      "Gene12",
      "Gene13",
      "Gene14",
      "Gene15",
      "Gene16",
      "Gene17",
      "Gene18",
      "Gene19",
      "Gene20"
    ]
  }
  ```

- **Ground Truth**: Marker lists specified in each JSON (e.g., PECAM1/CLDN5/KDR/FLT1/KLF2 for endothelial).
- **Tolerance**: Precision@20 must be ≥ 0.60 and Recall@20 must be ≥ 0.50.
- **Core Test**: Confirms differential expression rankings emphasise canonical markers.

### Canonical marker separation panels

- **Evals**: `*_canonical_expression_separation_v1.json`

  - [Endothelial](cell_typing/mouse_brain_endothelial_canonical_expression_separation_v1.json)
  - [GABAergic interneuron](cell_typing/mouse_brain_gabaergic_canonical_expression_separation_v1.json)
  - [Glutamatergic neuron](cell_typing/mouse_brain_glutamatergic_canonical_expression_separation_v1.json)
  - [Microglia](cell_typing/mouse_brain_microglia_canonical_expression_separation_v1.json)
  - [Neural progenitor](cell_typing/mouse_brain_neural_progenitor_canonical_expression_separation_v1.json)
  - [Oligodendrocyte](cell_typing/mouse_brain_oligodendrocyte_canonical_expression_separation_v1.json)

- **What & Why**: AUROC-based separation shows how cleanly canonical marker genes distinguish the target cluster from all others, strengthening biological interpretability.
- **Goal**: Compute AUROC for each canonical marker and report the per-gene scores plus the mean.
- **Input**: `latch://38438.account/AtlasxOmics/Pieper_154_brain_ArchRFull_10Core/annotated/D01900_NG05698_TBI_24hr.annotated.h5ad`
- **Expected Output**:

  ```json
  {
    "per_gene_stats": [
      {"gene": "GeneA", "auroc": 0.95},
      {"gene": "GeneB", "auroc": 0.92},
      {"gene": "GeneC", "auroc": 0.90},
      {"gene": "GeneD", "auroc": 0.88},
      {"gene": "GeneE", "auroc": 0.86},
      {"gene": "GeneF", "auroc": 0.84}
    ],
    "mean_auroc": 0.91
  }
  ```

- **Ground Truth**: Mean AUROC ≥ 0.85; at least 70 % of markers must have AUROC ≥ 0.80.
- **Tolerance**: Mean AUROC must be ≥ 0.85 and ≥ 70 % of markers must reach AUROC ≥ 0.80.
- **Core Test**: Demonstrates that canonical markers show strong expression separation in the identified cluster.

---

## Differential Accessibility (`atlasxomics/de`)
