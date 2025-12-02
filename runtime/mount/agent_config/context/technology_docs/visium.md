# Visium

0. Invoke the `redeem_package` tool to install the necessary Visium tools into the workspace.
  - `package_code`: `<PLACEHOLDER/>`
  - `package_version_id`: `<PLACEHOLDER/>`

1. Identify which Visium assay and chemistry the user is working with.  

2. Collect the appropriate Space Ranger output directory.  

3. Run `wf.__init__.visium_spaceranger_to_h5ad_workflow`, then load and visualize the resulting H5AD.  

4. Inspect the AnnData object to determine which preprocessing steps Space Ranger has already performed and skip accordingly.  

5. Choose the preprocessing backend for QC, normalization, PCA, integration, UMAP, clustering:

   **5.1. Assess dataset size**
   - Read `adata.n_obs` (spots) and `adata.n_vars` (genes)
   - Optionally check on-disk H5AD size

   **5.2. Use standard `scanpy` when**
   - `adata.n_obs ≤ 100,000`
   - H5AD size < ~10 GiB
   - Expected CPU preprocessing time is within tens of minutes
   - User does not need to re-run preprocessing many times

   **5.3. Use RAPIDS when**
   - `adata.n_obs > 100,000` or H5AD ≥ ~10 GiB  
   - CPU preprocessing is likely to be slow or already timing out
   - User expects to iterate on parameters frequently and needs faster re-runs  
   - Examples: multi-slide Visium integrations, Visium HD, large multi-sample projects

   **If RAPIDS is selected**
   - Run `wf.__init__.rapids_single-cell_preprocessing` on the H5AD (See `technology_docs/rapids.md`)
   - Continue downstream analysis using the RAPIDS-processed H5AD  

6. **Perform differential expression (1-vs-all):**
   - **Use RAPIDS DE** when `adata.n_obs > 100,000`  
   - **Use scanpy DE** when `adata.n_obs ≤ 100,000`  

7. Annotate each cluster by evaluating its top DE marker genes and assigning the best-matching cell type.

## Decision Tree: When to Use RAPIDS

START: Do you have a Visium HD dataset?
├─ YES, n_obs > 100K
│  ├─ Is it raw SpaceRanger output (no PCA/UMAP)?
│  │  └─ YES → RAPIDS full workflow
│  │
│  ├─ Is it already preprocessed (has PCA/UMAP/clustering)?
│  │  ├─ Need differential expression? 
│  │  │  └─ YES → RAPIDS with skip_qc, skip_normalization, skip_pca, skip_umap, skip_clustering=True
│  │  │  └─ NO → Manual annotation, visualization only
│  │  │
│  │  └─ Rerun with different parameters?
│  │     └─ YES → RAPIDS (iterate faster)
│  
└─ NO, n_obs ≤ 100K → Use standard scanpy

---

<about>

10x Genomics Visium is a suite of spatial transcriptomics kits, each optimized for different tissue types, input formats, and resolution needs. The main Visium product lines include:

- **Visium HD** — High-density barcoded array with ~2 µm pixels for near–single-cell resolution. Supports whole-transcriptome profiling from FFPE tissue using an enhanced chemistry and redesigned slide layout.

- **Visium for FFPE (v1 chemistry)** — Whole-transcriptome spatial profiling for FFPE tissues using probe-based capture. Provides ~55 µm spot resolution and is widely used for clinical and archival samples.

- **Visium for Fresh Frozen (v1 chemistry)** — Poly-A–based capture optimized for fresh-frozen tissues. Enables unbiased whole-transcriptome spatial analysis with standard Visium spot resolution.

- **Visium CytAssist Workflow (v2 chemistry)** — Updated Visium chemistry where tissue is mounted on standard glass slides and transferred to the capture slide using the CytAssist instrument. Improves tissue handling, compatibility, and standardization across labs.

- **Visium Gene Expression 3’ Kit (early Visium 3’/v1)** — The original Visium 3’ whole-transcriptome workflow for fresh-frozen samples, superseded by newer FFPE/v2 and HD kits but still used in legacy datasets.

These kits collectively span fresh-frozen, FFPE, and high-resolution use cases, giving researchers flexibility in tissue compatibility, experimental design, and spatial resolution.

</about>

---

<spaceranger_outputs>

# Understanding Space Ranger Outputs

Users often take raw blablabla from machine and run it thorugh 10x cloud blabla

## Overview of output structure  
The `spaceranger count` and `spaceranger aggr` pipelines run in a directory named after the `--id` argument (the "pipestance" directory). Output files appear in the `outs/` subdirectory within the pipestance directory.

The exact output files produced for a given pipestance depend on:

- the Space Ranger version used  
- which pipeline was used (`spaceranger count` or `spaceranger aggr`)  
- the parameters specified to the pipeline  
- whether or not the CytAssist instrument was used  
- whether brightfield or fluorescence images were provided  
- which Visium assay was used  

### Visium HD and Visium HD 3′ Gene Expression  
When `spaceranger count` is run on Visium HD or Visium HD 3′ data, the outputs have a hierarchical structure to accommodate cell segmentation and various binning levels.

| File or Directory Name | Description |
|------------------------|-------------|
| `barcode_mappings.parquet` | This file efficiently stores spatial mapping information, essentially functioning as a CSV that tracks relationships between barcodes (squares), nuclei, cells, and bins within your Visium HD data. |
| `binned_outputs` | By default, this directory has three subdirectories: `square_002um`, `square_008um`, and `square_016um`. Each directory contains `filtered_feature_bc_matrix`, `raw_feature_bc_matrix`, `spatial`, `filtered_feature_bc_matrix.h5`, and `raw_feature_bc_matrix.h5`. The `analysis` directory is only provided at 8 and 16 µm bin size. The `cloupe.cloupe` is only provided at 8 µm bin size. The `raw_probe_bc_matrix.h5` is only provided at 2 µm resolution. |
| `cloupe_008um.cloupe` | Symlink to the .cloupe file at 8 µm bin size |
| `cloupe_cell.cloupe` | Symlink to the .cloupe file with cell segmentation |
| `feature_slice.h5` | A new file type, specific to Visium HD, to support efficient fetching of 2 µm resolution image slices for a single gene or multiple genes. |

Other files:  
- `metrics_summary.csv`: Run summary metrics in CSV format  
- `molecule_info.h5`: Contains per-molecule information for all molecules that contain a valid barcode, valid UMI, and were assigned with high confidence to a gene barcode or bin.  
- `probe_set.csv`: Copy of the input probe set reference CSV file.  
- `segmented_outputs`: Folder containing segmented outputs. Contains `analysis`, `cell_segmentations.geojson`, `cloupe.cloupe`, `filtered_feature_cell_matrix`, `filtered_feature_cell_matrix.h5`, `graphclust_annotated_cell_segmentations.geojson`, `graphclust_annotated_nucleus_segmentations.geojson`, `nucleus_segmentations.geojson`, `raw_feature_cell_matrix`, `raw_feature_cell_matrix.h5`, and `spatial`.  
- `spatial`: Folder containing outputs that capture the spatiality of the data.  
- `web_summary.html`: Run summary metrics and plots in HTML format.  

**Example SpaceRanger Output Folder for Visium HD and Visium HD 3′ Gene Expression**

visium_hd_data
|-- binned_outputs
|   |-- square_008um
|   |   |-- analysis
|   |   |   |-- clustering
|   |   |   |   |-- gene_expression_graphclust
|   |   |   |   |   `-- clusters.csv
|   |   |   |-- diffexp
|   |   |   |   |-- gene_expression_graphclust
|   |   |   |   |   `-- differential_expression.csv
|   |   |   |-- pca
|   |   |   |   `-- gene_expression_10_components
|   |   |   |       |-- components.csv
|   |   |   |       |-- dispersion.csv
|   |   |   |       |-- features_selected.csv
|   |   |   |       |-- projection.csv
|   |   |   |       `-- variance.csv
|   |   |   `-- umap
|   |   |       `-- gene_expression_2_components
|   |   |           `-- projection.csv
|   |   |-- cloupe.cloupe
|   |   |-- filtered_feature_bc_matrix
|   |   |   |-- barcodes.tsv.gz
|   |   |   |-- features.tsv.gz
|   |   |   `-- matrix.mtx.gz
|   |   |-- filtered_feature_bc_matrix.h5
|   |   |-- raw_feature_bc_matrix
|   |   |   |-- barcodes.tsv.gz
|   |   |   |-- features.tsv.gz
|   |   |   `-- matrix.mtx.gz
|   |   |-- raw_feature_bc_matrix.h5
|   |   `-- spatial
|   |       |-- aligned_fiducials.jpg
|   |       |-- aligned_tissue_image.jpg
|   |       |-- cytassist_image.tiff
|   |       |-- detected_tissue_image.jpg
|   |       |-- scalefactors_json.json
|   |       |-- tissue_hires_image.png
|   |       |-- tissue_lowres_image.png
|   |       `-- tissue_positions.parquet
`-- spatial
    |-- aligned_fiducials.jpg
    |-- aligned_tissue_image.jpg
    |-- cytassist_image.tiff
    |-- detected_tissue_image.jpg
    |-- tissue_hires_image.png
    `-- tissue_lowres_image.png

### Visium v1/v2 Spatial Gene Expression  
When running `spaceranger count` on Visium v1/v2 Spatial Gene Expression libraries, the following files can be found within the `outs/` subfolder:

- `web_summary.html`: Run summary metrics and plots in HTML format  
- `cloupe.cloupe`: Loupe Browser visualization and analysis file  
- `spatial/`: Folder containing outputs that capture the spatiality of the data.  
- `analysis/`: Folder containing secondary analysis data including graph-based clustering and K-means clustering (K = 2–10); differential gene expression between clusters; PCA, t-SNE, and UMAP dimensionality reduction.  
- `metrics_summary.csv`: Run summary metrics in CSV format  
- `probe_set.csv`: Copy of the input probe set reference CSV file. Present for Visium FFPE and CytAssist workflow  
- `possorted_genome_bam.bam`: Indexed BAM file containing position‐sorted reads aligned to the genome and transcriptome, annotated with barcode information  
- `possorted_genome_bam.bam.bai`: Index for `possorted_genome_bam.bam`. In cases where the reference transcriptome is generated from a genome with very long chromosomes (> 512 Mbp), Space Ranger v2.0+ generates a `possorted_genome_bam.bam.csi` index file instead.  
- `filtered_feature_bc_matrix/`: Contains only tissue-associated barcodes in MEX format. Each element of the matrix is the number of UMIs associated with a feature (row) and a barcode (column). This file can be input into third-party packages and allows users to wrangle the barcode-feature matrix (e.g., to filter outlier spots, run dimensionality reduction).  
- `filtered_feature_bc_matrix.h5`: Same information as `filtered_feature_bc_matrix/` but in HDF5 format.  
- `raw_feature_bc_matrices/`: Contains all detected barcodes in MEX format. Each element of the matrix is the number of UMIs associated with a feature (row) and a barcode (column).  
- `raw_feature_bc_matrix.h5`: Same information as `raw_feature_bc_matrices/` but in HDF5 format.  
- `raw_probe_bc_matrix.h5`: Contains UMI counts of each probe for all detected barcodes in HDF5 format. Only produced when running pipelines for probe-based assays.  
- `molecule_info.h5`: Contains per-molecule information for all molecules that contain a valid barcode, valid UMI, and were assigned with high confidence to a gene or protein barcode. This file is required for additional analysis `spaceranger` pipelines including `aggr`, `targeted-compare` and `targeted-depth`.  
- Space Ranger v2.1 introduced support for Gene and Protein Expression (v2 slides only), and there are a few variations from the list above.

</spaceranger_outputs>

---

<spaceranger_to_anndata>

## SpaceRanger to H5AD Workflow

**Currently only supports Visium HD and Visium 3' Space Ranger outputs**
- Launch `wf.__init__.visium_spaceranger_to_h5ad_workflow` and wait for execution outputs. 

```python
from lplots.widgets.workflow import w_workflow

w = w_workflow(
    label="Convert SpaceRanger Outputs to H5AD",
    wf_name="wf.__init__.visium_spaceranger_to_h5ad_workflow",
    version=None, # Uses latest version
    params={
      "input_file": LatchDir("latch://38438.account/visium/colon"),
      "bin_name": "square_008um",
      "output_directory": LatchOutputDir("latch://38438.account/spaceranger_outs/test_colon"),
      "run_name": "colon_test"
    }
)

execution = w.value

if execution is not None:
    res = await execution.wait()
    workflow_outputs = list(res.output.values())

    if len(workflow_outputs) == 0:
        w_text_output(content="No workflow outputs")
        exit()
    
    output_dir = workflow_outputs[0] 
    
    browser = w_ldata_browser(
        label="Browse Workflow Outputs", 
        dir=output_dir
    )

    h5ad_file = None
    
    for item in output_dir.iterdir():
        item_name = item.name()
        if item_name and item_name.endswith('.h5ad'):
            h5ad_file = item

    # Download the H5AD file + visualize using w_h5

```

</spaceranger_to_anndata>
