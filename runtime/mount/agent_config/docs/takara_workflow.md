## Takara Seeker and Trekker Analysis Workflow

This is the **authoritative step-by-step pipeline** for Takara Seeker and Trekker experiments.
Follow steps in order. Apply Seeker-specific preprocessing when required:

1. **Data Loading** - load data using **Scanpy**
2. **Experiment Setup** - Ask users to confirm whether their experiment is Seeker or Trekker, and if their H5AD contains one or multiple samples.
3. **Background Removal** - *Seeker only*
4. **Quality Control & Filtering** - use **Scanpy** to remove low-quality cells/spots
5. **Normalization** - e.g., total-count scaling, log1p
6. **Feature Selection** - identify highly variable genes (HVGs)
7. **Dimensionality Reduction** - compute **PCA**, then **UMAP** embeddings
8. **Clustering** - apply **Leiden** clustering
9. **Differential Gene Expression (DGE)** - find cluster marker genes
10. **Cell Type Annotation** - assign biological meaning to clusters

Some pointers on steps.

### Experiment Setup

- You must ask users to confirm whether their experiment is Seeker 3x3, Seeker 10x10, or Trekker.
- The user's answer will directly influence the code written for background removal.


### Background Removal

*Goal: Remove off-tissue beads (from Curio Seeker data ONLY) to create an on-tissue mask for downstream analysis.*

Three-step filter (applied sequentially):

1/ UMI threshold: keep beads with log10(UMI) ≥ t.
- Pick t at the local minimum between the two modes in the log10(UMI) histogram.
- Typical range: 1–2; expose as widget min_log10_UMI.

2/ Neighborhood density A: square window m×m μm (default m=40).
- Keep beads with count ≥ p (default p=5).

3/ Neighborhood density B: square window n×n μm (default n=100).
- Keep beads with count ≥ q (default q=10).


UI & Plots to generate:

- Histogram of log10(UMI) with adjustable t and vertical threshold line.
- Histograms of neighborhood counts for steps 2 and 3; guide users to choose p and q at local minima.
- Spatial scatter/overlay showing kept vs removed beads after each step.
- Expose widgets for t, m, n, p, q; default to t=1.6, m=40, p=5, n=100, q=10.
