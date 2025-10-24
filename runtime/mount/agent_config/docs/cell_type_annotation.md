# Tangram-Based Spatial Cell Type Annotation

## Purpose

This workflow provides **spatial cell type annotation** using **Tangram** with pre-aggregated **Tabula Sapiens** reference data. Tangram maps reference expression profiles onto spatial data to predict cell type identities and project gene expression patterns.

## When to Use

- **Human spatial data only** (CosMX, Visium, Slide-seq, Xenium, MERFISH, etc.)
- When you need **continuous cell type probabilities** rather than discrete labels
- When you want to **project full expression profiles** onto spatial coordinates
- For **validation** of marker gene-based annotations

## Advantages & Limitations

**Advantages:**

- Pre-aggregated clusters = **faster training**, less memory
- Comprehensive human reference from Tabula Sapiens
- Provides **cell type probabilities** for each spatial location
- Can project **entire transcriptome** onto spatial data
- Validated method published in [Nature Methods](https://www.nature.com/articles/s41592-021-01264-7)

**Limitations:**

- **Human data only** (Tabula Sapiens is human-specific)
- **Requires raw count data** — both reference and spatial must be in raw counts (not normalized)
- Requires reasonable **gene overlap** between reference and spatial data
- Training can take time (10-30 minutes depending on data size)
- Works best when reference and spatial data are from similar tissues

---

## Data Overview

### Tabula Sapiens Reference Profiles

Pre-aggregated reference dataset in **raw counts**:

**`latch://{account_id}/tabula_sapiens_avg_raw.csv.gz`**

- Raw count data averaged across cell types
- **Use raw counts for Tangram** — Tangram works best when reference and spatial data are in the same format
- Spatial data should also be in raw counts (not normalized)

**Structure:**

- **Index**: Multi-index tuples of `(tissue, cell_type)`
- **Columns**: Gene symbols (human)
- **Values**: Average expression for that gene in that cell type/tissue combination

**Example index:**

```python
('Blood', 'B cell')
('Blood', 'T cell')
('Brain', 'astrocyte')
('Brain', 'excitatory neuron')
('Liver', 'hepatocyte')
...
```

---

## Workflow

### CRITICAL: Widget Selection Pattern for Tangram

**Tangram training is time-consuming (10-30 minutes).** Therefore:

1. ✅ **Create widgets to collect parameters** (device, epochs, tissue selection)
2. ✅ **WAIT for user to make selections** — do NOT proceed automatically
3. ✅ **Split into separate cells**:
   - **Cell A**: Create widgets and display them
   - **Cell B**: Run Tangram training using selected values
4. ❌ **Do NOT run Tangram immediately** after creating widgets in the same cell

**Why this matters:**

- Tangram takes 10-30 minutes to train
- Users need to confirm parameters before committing to training
- If parameters are wrong, users waste time re-running
- This is different from marker gene annotation (which is fast and can run automatically)

**Correct pattern:**

```python
# Cell 1: Create widgets and wait
device = w_select(label="Device", options=["cpu", "cuda:0"])
epochs = w_number_input(label="Epochs", default=500)
print("Please select parameters above, then run the next cell.")
```

```python
# Cell 2: Only run after user confirms (separate cell)
ad_map = tg.map_cells_to_space(
    adata_sc=ad_ref,
    adata_sp=ad_sp,
    device=device.value,
    num_epochs=epochs.value
)
```

**Incorrect pattern:**

```python
# ❌ Don't do this - runs immediately without user confirmation
device = w_select(label="Device", options=["cpu", "cuda:0"])
epochs = w_number_input(label="Epochs", default=500)

# Tangram starts immediately - user has no chance to adjust!
ad_map = tg.map_cells_to_space(...)  # Takes 30 minutes!
```

---

### Step 1: Load Tabula Sapiens Reference

Download and convert the reference profiles to AnnData format for Tangram compatibility.

```python
import pandas as pd
import numpy as np
import anndata as ad
from latch.ldata.path import LPath
from pathlib import Path

# Download reference (raw counts)
ref_lpath = LPath("latch://{account_id}/tabula_sapiens_avg_raw.csv.gz")
local_ref_path = Path(f"{ref_lpath.node_id()}.csv.gz")
ref_lpath.download(local_ref_path, cache=True)

# Load reference DataFrame
ref_df = pd.read_csv(local_ref_path, sep='\t', index_col=0)

# Parse multi-index: index strings are "(tissue, cell_type)" tuples
# Convert string representation to actual tuples
import ast
tissues = []
cell_types = []

for idx in ref_df.index:
    tissue, cell_type = ast.literal_eval(idx)
    tissues.append(tissue)
    cell_types.append(cell_type)

# Create AnnData from reference
# Each "cell" is actually an aggregated cell type profile
ad_ref = ad.AnnData(X=ref_df.values, var=pd.DataFrame(index=ref_df.columns))
ad_ref.obs['tissue'] = tissues
ad_ref.obs['cell_type'] = cell_types

# Store combined annotation for Tangram cluster mode
ad_ref.obs['tissue_cell_type'] = [f"{t}_{ct}" for t, ct in zip(tissues, cell_types)]

print(f"Reference: {ad_ref.n_obs} cell types × {ad_ref.n_vars} genes")

# CRITICAL: Introspect available tissues from the data
available_tissues = sorted(pd.unique(ad_ref.obs['tissue']))
print(f"Available tissues in reference: {available_tissues}")
print(f"Number of tissues: {len(available_tissues)}")
```

### Step 2: Prepare Spatial Data

Load spatial data and ensure it's in **raw counts** to match the reference.

**CRITICAL:** Tangram requires that both reference and spatial data are in the same format. Since the reference is raw counts, spatial data must also be raw counts.

```python
import scanpy as sc
import numpy as np

# Load spatial AnnData (example: from h5ad file)
ad_sp = sc.read_h5ad("spatial_data.h5ad")

# Ensure spatial coordinates exist
if 'spatial' not in ad_sp.obsm:
    raise ValueError("Spatial data must have coordinates in adata.obsm['spatial']")

# CRITICAL: Check if spatial data is in raw counts
# Raw counts should be integers or close to integers
is_raw_counts = np.all(ad_sp.X.data == ad_sp.X.data.astype(int))

if not is_raw_counts:
    print("WARNING: Spatial data does not appear to be raw counts!")
    print("Tangram works best when both reference and spatial data are in raw counts.")
    print("If you have raw counts available, use adata.raw or reload from raw data.")

    # If raw counts are stored in adata.raw
    if ad_sp.raw is not None:
        print("Using raw counts from adata.raw...")
        ad_sp = ad_sp.raw.to_adata()
    else:
        print("No raw counts found. Results may be suboptimal.")
        print("Consider reloading data without normalization.")

# Check gene overlap
common_genes = list(set(ad_ref.var_names) & set(ad_sp.var_names))
print(f"Gene overlap: {len(common_genes)} / {ad_ref.n_vars} reference genes")

if len(common_genes) < 100:
    print(f"WARNING: Low gene overlap ({len(common_genes)} genes). Results may be unreliable.")
    print("Consider using marker gene annotation instead.")
```

### Step 3: Filter to Common Genes and Prepare for Tangram

Use Tangram's preprocessing function to align the datasets.

```python
import tangram as tg

# Option 1: Use all common genes (recommended for pre-aggregated reference)
tg.pp_adatas(ad_ref, ad_sp, genes=None)

# Option 2: Use specific training genes (if you have marker genes)
# training_genes = ['CD3D', 'CD8A', 'MS4A1', 'CD14', ...]  # Your marker genes
# tg.pp_adatas(ad_ref, ad_sp, genes=training_genes)

# Verify alignment
assert ad_ref.uns['training_genes'] == ad_sp.uns['training_genes']
print(f"Training genes: {len(ad_ref.uns['training_genes'])}")
```

### Step 4: Map Cell Types to Space (Cluster-Level Tangram)

Run Tangram mapping using the pre-aggregated reference in cluster mode.

**This is the key step**: Tangram learns to map reference cell types onto spatial coordinates.

**CRITICAL:** This step should be split into TWO cells:

1. **Cell 1**: Create widgets and let user make selections
2. **Cell 2**: Run Tangram training with selected parameters

#### Cell 1: Collect User Parameters (Create Widgets)

```python
from lplots.widgets.select import w_select
from lplots.widgets.number import w_number_input

# Widget for choosing device
device_select = w_select(
    label="Device for training (CPU is slower but uses less memory)",
    options=["cpu", "cuda:0"],
    default="cpu"
)

# Widget for number of training epochs
epochs_input = w_number_input(
    label="Number of training epochs (500-1000 recommended)",
    default=500,
    min_value=100,
    max_value=2000,
    step=100
)

# Optional: Widget for tissue filtering
from lplots.widgets.select import w_select
import pandas as pd

available_tissues = sorted(pd.unique(ad_ref.obs['tissue']))
tissue_filter = w_select(
    label="Filter reference to specific tissue (optional)",
    options=["Use all tissues"] + available_tissues,
    default="Use all tissues"
)

print("Please select your parameters above, then proceed to the next cell to start training.")
print(f"Note: Tangram training will take approximately 10-30 minutes on CPU.")
```

**STOP HERE** — Wait for user to make selections in the widgets above before proceeding.

#### Cell 2: Run Tangram Training (After User Selections)

**Only run this cell after the user has made their selections above.**

```python
import tangram as tg

# Filter reference by tissue if specified
if tissue_filter.value != "Use all tissues":
    ad_ref_filtered = ad_ref[ad_ref.obs['tissue'] == tissue_filter.value].copy()
    print(f"Using {ad_ref_filtered.n_obs} cell types from {tissue_filter.value}")
else:
    ad_ref_filtered = ad_ref
    print(f"Using all {ad_ref_filtered.n_obs} cell types from all tissues")

# Map cells to space using cluster mode
# This uses the pre-aggregated cell type profiles
print(f"Starting Tangram training with {epochs_input.value} epochs on {device_select.value}...")
print("This will take 10-30 minutes. Please wait...")

ad_map = tg.map_cells_to_space(
    adata_sc=ad_ref_filtered,
    adata_sp=ad_sp,
    mode='clusters',  # Use cluster mode for pre-aggregated data
    cluster_label='cell_type',  # Field in ad_ref.obs containing cell type labels
    device=device_select.value,
    num_epochs=epochs_input.value
)

print("✓ Mapping complete!")
print(f"Final mapping score: {ad_map.uns.get('training_history', {}).get('train_score', 'N/A')}")
```

**Training notes:**

- Mapping typically takes 10-30 minutes on CPU
- Use GPU (`cuda:0`) for faster training if available
- Higher scores indicate better alignment
- Score should plateau — if it's still increasing, consider more epochs

### Step 5: Project Cell Annotations (KEY OUTPUT)

**This is the PRIMARY output:** Transfer cell type labels from reference to spatial data.

```python
# Project cell type annotations onto spatial coordinates
tg.project_cell_annotations(ad_map, ad_sp, annotation='cell_type')

print("Cell type annotations projected to spatial data!")
print(f"Annotations stored in: ad_sp.obsm['tangram_ct_pred']")

# The result is a cells × cell_types probability matrix
# Each row sums to 1 (probability distribution over cell types)
cell_type_probs = ad_sp.obsm['tangram_ct_pred']
print(f"Probability matrix shape: {cell_type_probs.shape}")

# Get most likely cell type for each spatial location
predicted_cell_types = cell_type_probs.columns[cell_type_probs.values.argmax(axis=1)]
ad_sp.obs['tangram_cell_type'] = predicted_cell_types

# Get prediction confidence (max probability)
ad_sp.obs['tangram_confidence'] = cell_type_probs.values.max(axis=1)
```

### Step 6: Visualize Mapped Cell Types

Display cell type predictions on spatial coordinates using `w_h5`.

```python
from lplots.widgets.h5 import w_h5
from lplots.widgets.table import w_table
import pandas as pd

# Visualize with w_h5
viewer = w_h5(ann_data=ad_sp)

# Summary table: cell type distribution
cell_type_counts = ad_sp.obs['tangram_cell_type'].value_counts().reset_index()
cell_type_counts.columns = ['Cell Type', 'Count']
cell_type_counts['Percentage'] = 100 * cell_type_counts['Count'] / cell_type_counts['Count'].sum()

w_table(source=cell_type_counts, label="Tangram Cell Type Distribution")

# Confidence summary
confidence_stats = ad_sp.obs['tangram_confidence'].describe().to_frame('Confidence')
w_table(source=confidence_stats, label="Prediction Confidence Statistics")
```

### Step 7: Project Gene Expression (Optional)

Optionally, project the full reference transcriptome onto spatial coordinates.

**Use case:** Visualize expected gene expression patterns based on predicted cell types.

```python
# Project genes from reference onto spatial data
ad_ge = tg.project_genes(
    adata_map=ad_map,
    adata_sc=ad_ref,
    cluster_label='cell_type'
)

# ad_ge is a new AnnData with:
# - Same spatial coordinates as ad_sp
# - Gene expression projected from reference

# Example: visualize a projected gene
gene_to_plot = 'CD3D'  # T cell marker
if gene_to_plot in ad_ge.var_names:
    ad_ge.obs[f'{gene_to_plot}_projected'] = ad_ge[:, gene_to_plot].X.toarray().flatten()
    viewer_ge = w_h5(ann_data=ad_ge)
```

---

## Integration with Existing Workflows

### Relationship to Marker Gene Annotation

Tangram and marker gene annotation are **complementary approaches**:

| Feature             | Marker Gene Annotation        | Tangram Mapping                 |
| ------------------- | ----------------------------- | ------------------------------- |
| **Data**            | Any organism                  | Human only                      |
| **Speed**           | Fast (< 1 min)                | Slower (10-30 min)              |
| **Output**          | Discrete labels               | Continuous probabilities        |
| **Gene projection** | No                            | Yes                             |
| **Method**          | Database lookup               | ML alignment                    |
| **Best for**        | Quick annotation, any species | Detailed human spatial analysis |

### Combined Approach: Use Both for Validation

**Recommended workflow for human data:**

1. **Run marker gene annotation first** (fast, discrete labels)
2. **Run Tangram mapping second** (detailed probabilities)
3. **Compare results** to validate and refine annotations

```python
# Compare Tangram vs marker gene predictions
if 'predicted_cell_type' in ad_sp.obs and 'tangram_cell_type' in ad_sp.obs:
    comparison = pd.DataFrame({
        'Marker Gene Prediction': ad_sp.obs['predicted_cell_type'],
        'Tangram Prediction': ad_sp.obs['tangram_cell_type'],
        'Tangram Confidence': ad_sp.obs['tangram_confidence']
    })

    # Find agreements
    comparison['Agreement'] = (
        comparison['Marker Gene Prediction'] == comparison['Tangram Prediction']
    )

    agreement_rate = 100 * comparison['Agreement'].mean()
    print(f"Agreement rate: {agreement_rate:.1f}%")

    # Display disagreements (interesting cases)
    disagreements = comparison[~comparison['Agreement']]
    if len(disagreements) > 0:
        w_table(
            source=disagreements.head(20),
            label="Top 20 Prediction Disagreements (investigate these)"
        )
```

---

## Troubleshooting

### Low Gene Overlap

**Problem:** Fewer than 100-200 common genes between reference and spatial data.

**Solutions:**

- Check gene naming (uppercase vs lowercase)
- Convert gene symbols to common format (e.g., all uppercase)
- If overlap is < 100 genes, consider using marker gene annotation instead
- For mouse data, use marker gene annotation (Tangram reference is human-only)

```python
# Check and fix gene naming
ad_sp.var_names = ad_sp.var_names.str.upper()  # Convert to uppercase
ad_ref.var_names = ad_ref.var_names.str.upper()

# Re-check overlap
common_genes = list(set(ad_ref.var_names) & set(ad_sp.var_names))
print(f"Gene overlap after fixing: {len(common_genes)}")
```

### Low Mapping Score / Poor Convergence

**Problem:** Training score plateaus at < 0.5 or keeps decreasing.

**Solutions:**

- Increase training epochs
- Check that spatial data has reasonable gene expression (not all zeros)
- Verify spatial and reference data are from compatible tissues
- Try using marker genes instead of all genes for training

```python
# Check spatial data quality
print(f"Non-zero genes per cell: {(ad_sp.X.toarray() > 0).sum(axis=1).mean():.1f}")
print(f"Cells with < 50 genes: {((ad_sp.X.toarray() > 0).sum(axis=1) < 50).sum()}")

# Try with marker genes only
from lplots.widgets.select import w_multiselect

# Let user select tissue-specific markers
# ...then run tg.pp_adatas with those markers
```

### Memory Issues

**Problem:** Out of memory during training.

**Solutions:**

- Use CPU instead of GPU (lower memory but slower)
- Reduce number of spatial cells (subsample)
- Use fewer training genes
- Close other applications

```python
# Subsample spatial data if too large
if ad_sp.n_obs > 50000:
    print(f"Subsampling from {ad_sp.n_obs} to 50000 cells...")
    sc.pp.subsample(ad_sp, n_obs=50000, random_state=42)
```

### Tissue Mismatch

**Problem:** Reference tissue doesn't match spatial data tissue, or you want to use only a specific tissue subset.

**Solutions:**

- **Always introspect available tissues first** — extract from data, don't hard-code tissue names
- Filter reference to relevant tissue only
- Let user select from available options using widget
- Validate results carefully with marker genes

```python
# CRITICAL: Always introspect available tissues from the data first
from lplots.widgets.select import w_select
import pandas as pd

# Extract available tissues from the reference data
available_tissues = sorted(pd.unique(ad_ref.obs['tissue']))
print(f"Available tissues in reference: {available_tissues}")

# Present tissue options to user for selection
tissue_selector = w_select(
    label="Select reference tissue to use for mapping",
    options=available_tissues,
    default=available_tissues[0] if available_tissues else None
)

# Filter reference to selected tissue
if tissue_selector.value:
    ad_ref_filtered = ad_ref[ad_ref.obs['tissue'] == tissue_selector.value].copy()
    print(f"Filtered reference: {ad_ref_filtered.n_obs} cell types from {tissue_selector.value}")

    # Now use ad_ref_filtered for mapping
else:
    print("No tissue selected, using full reference")
    ad_ref_filtered = ad_ref
```

**Important:** Never hard-code tissue names like "Brain", "Liver", etc. Always extract from the data using `pd.unique(ad_ref.obs['tissue'])` to see what's actually available.

---

## Best Practices

1. **CRITICAL: Split widget creation from Tangram training** — create widgets in one cell, run training in the next cell after user selections
2. **Wait for user parameter selection** — Tangram takes 10-30 minutes, users must confirm parameters first
3. **Use raw counts for both reference and spatial data** — CRITICAL for Tangram alignment
4. **Verify spatial data has raw counts and use them** — check before mapping
5. **Always introspect available tissues** — extract from data using `pd.unique(ad_ref.obs['tissue'])`, never hard-code tissue names
6. **Check gene overlap first** — need at least 100-200 common genes
7. **Match tissues** — filter reference to relevant tissue when possible using extracted tissue list
8. **Monitor training** — check that score plateaus and doesn't keep decreasing
9. **Validate results** — compare with marker gene annotations when available
10. **Check confidence** — low confidence predictions may need manual review
11. **Use cluster mode** — pre-aggregated reference is faster and uses less memory
12. **Cache downloads** — download reference once per session for efficiency
13. **Start with CPU** — use CPU first, switch to GPU only if needed
14. **Visualize extensively** — use `w_h5` to inspect spatial patterns

---

## Summary

Tangram-based spatial cell type annotation provides:

- **Probabilistic cell type assignments** for human spatial data
- **Gene expression projection** from comprehensive reference
- **Validation** of marker gene-based predictions
- **Granular cell type resolution** from Tabula Sapiens

**Workflow summary:**

1. Load Tabula Sapiens reference (**raw counts**)
2. Load and prepare spatial data (**ensure raw counts**)
3. Map reference cell types to space using cluster-level Tangram
4. **Project cell annotations** (key output: cell type probabilities)
5. Optionally project gene expression
6. Visualize and validate results

**When to use:**

- Human spatial data with good gene coverage
- Need for continuous cell type probabilities
- Want to project reference expression patterns
- Validation of discrete annotations

**Complementary to marker gene annotation** — use both for comprehensive analysis!
