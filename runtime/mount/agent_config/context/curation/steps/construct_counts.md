<goal>
Parse downloaded supplementary files into a standardized AnnData object with Ensembl gene IDs, sample identifiers, and raw counts.
</goal>

<context_files>
Read from download step:
- `/opt/latch/.../curation/tmp/{accession}/study_metadata.txt`
- `/opt/latch/.../curation/tmp/{accession}/paper_text.txt`
- `/opt/latch/.../curation/tmp/{accession}/downloaded_files.txt`

Write after construction:
- `/opt/latch/.../curation/tmp/{accession}/constructed_h5ad.txt` (path to output file)
</context_files>

<method>
**Setup:**
```python
import sys
from pathlib import Path

sys.path.insert(0, "/opt/latch/plots-faas/runtime/mount/agent_config/context/curation/lib")

CONTEXT_DIR = Path("/opt/latch/plots-faas/runtime/mount/agent_config/context/curation/tmp/{accession}")
```

**1. Inspect downloaded files:**
- Read `downloaded_files.txt` to get file paths
- Examine file types (.csv, .tsv, .h5ad, .mtx, etc.)
- Use paper_text and study_metadata for context on file contents

**2. Parse files into AnnData:**

For existing H5AD:
```python
import anndata as ad
adata = ad.read_h5ad(input_path)
```

For CSV/TSV count matrices:
```python
import pandas as pd
from anndata import AnnData

counts_df = pd.read_csv(counts_path, index_col=0)
adata = AnnData(X=counts_df.values, obs=pd.DataFrame(index=counts_df.index), var=pd.DataFrame(index=counts_df.columns))
```

For 10X MTX format:
```python
import scanpy as sc
adata = sc.read_10x_mtx(mtx_dir)
```

**3. Standardize gene IDs:**
```python
from curate import convert_and_swap_symbol_index, is_valid_ensembl

# Check if already Ensembl
if not all(is_valid_ensembl(g) for g in adata.var_names[:100]):
    adata = convert_and_swap_symbol_index(adata)
```

**4. Add sample identifier:**
```python
# Infer from metadata or filename
adata.obs["latch_sample_id"] = sample_name

# Prefix obs_names with sample to ensure uniqueness across samples
adata.obs_names = [f"{sample_name}_{bc}" for bc in adata.obs_names]
```

**5. Prefix author metadata:**
```python
# Rename existing obs columns
for col in adata.obs.columns:
    if col != "latch_sample_id" and not col.startswith("author_"):
        adata.obs.rename(columns={col: f"author_{col}"}, inplace=True)
```

**6. Merge multiple samples:**
```python
import anndata as ad
from curate import reindex_and_fill_list

# If gene sets differ, reindex to union
adatas = reindex_and_fill_list(adatas)
combined = ad.concat(adatas, join="outer")
```

**7. Validate:**
```python
from curate import validate_counts_object, format_validation_report, all_checks_passed

validation_log = validate_counts_object(adata)
print(format_validation_report(validation_log))

if not all_checks_passed(validation_log):
    # Fix issues and re-validate
    pass
```

**8. Save and record path:**
```python
output_path = Path("./output.h5ad")
adata.write_h5ad(output_path)
(CONTEXT_DIR / "constructed_h5ad.txt").write_text(str(output_path.resolve()))
```
</method>

<library>
curate.parsing
curate.validation
</library>

<validation_criteria>
The constructed AnnData must satisfy:
- `obs["latch_sample_id"]` exists
- `var.index` contains Ensembl IDs (pattern: `ENS[A-Z0-9]{0,5}G\d{11}`)
- `var["gene_symbols"]` exists with unique values
- `obs.index` is unique
- `var.index` is unique
- Counts are non-negative integers
- Matrix has non-zero values
- Any `author_*` columns are not all NaN
</validation_criteria>

<self_eval_criteria>
- All validation checks pass (use `all_checks_passed()`)
- Cell count roughly matches paper expectation (check paper_text)
- `constructed_h5ad.txt` written with valid path
- No data loss (didn't subset/sample/drop cells)
</self_eval_criteria>
