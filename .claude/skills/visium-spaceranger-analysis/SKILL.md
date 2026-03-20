---
name: visium-spaceranger-analysis
description: >
  Use this skill when the user is analyzing 10x Genomics Visium data, including
  Visium HD, CytAssist, FFPE, or Space Ranger outputs such as spatial/ folders,
  tissue_positions.csv, tissue_positions.parquet, binned_outputs, or
  filtered_feature_bc_matrix.h5. Use it for Space Ranger intake, H5AD conversion,
  preprocessing, RAPIDS routing for large datasets, differential expression,
  and cluster annotation.
---

# Visium Space Ranger Analysis

Use this skill for Visium platform detection, Space Ranger intake, preprocessing
routing, and downstream analysis order.

## Detect this platform when

- the user mentions `Visium`, `Visium HD`, `CytAssist`, or `Space Ranger`
- a directory contains `spatial/`
- files include `tissue_positions.csv` or `tissue_positions.parquet`
- files include `filtered_feature_bc_matrix.h5`
- a directory contains `binned_outputs`

## Latch setup

If the environment supports `redeem_package`, install the Visium package before analysis:

- `package_code`: `3e2e4bdf634b943cfd2571f22e5aaaae29bdce5206ad89c71090ee19d96df606`
- `package_version_id`: `402`

If `redeem_package` is unavailable, continue with the repo-local references and the current environment.

## Workflow overview

1. Identify the Visium assay and chemistry — see [Visium reference](references/visium_reference.md)
2. Collect the correct Space Ranger output directory — see [Visium reference](references/visium_reference.md)
3. Convert Space Ranger output to H5AD — see [Visium reference](references/visium_reference.md)
4. Inspect the AnnData object and skip preprocessing steps that were already completed
5. Choose Scanpy vs RAPIDS based on dataset size — see [RAPIDS preprocessing](references/rapids_preprocessing.md)
6. Run differential expression
7. Annotate clusters and perform a biological plausibility check

## Large dataset routing

- Use standard Scanpy when the dataset is modest and CPU preprocessing is practical.
- Use RAPIDS when the dataset is large, especially for Visium HD or multi-sample integrations.

## Latch-specific execution

If `latch-workflows`, `latch-plots-ui`, or `latch-data-access` are available, prefer them for:

- `w_workflow`
- `w_h5`, `w_plot`, `w_table`, and `w_text_output`
- Latch Data path selection and conversion to workflow inputs

If those sibling skills are not available, use the local reference docs directly.
