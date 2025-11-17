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

---

### Next Steps
- Refer to <atlasx_analysis_overview> to see the guide for other types of analyses. 
