## Differential Analysis Workflow

`compare_workflow` compares differences in genes, peaks, and motifs between user-defined cluster/condition groupings. Outputs include:

**For genes**:
- `volcano_plot.pdf`: Log2FC vs -log10(p-value) visualization
- `all_genes.csv`: Complete differential expression results from ArchR::getMarkerFeatures
- `marker_genes.csv`: Filtered significant genes with Log2FC, FDR, and enrichment scores

**For peaks**:
- `MA_plot.pdf`: Average signal vs Log2FC visualization
- `all_peaks.csv`: All peak accessibility differences
- `marker_peaks.csv`: Filtered significant peaks with genomic coordinates

**For motifs**:
- `[up/down]Regulated_motifs.csv`: TF motifs ranked by significance
- `[up/down]_enrichment_plot.pdf`: Scatter plots ranked by -log10(FDR)
- `all_motifs.csv`: Complete motif enrichment differences
- `marker_motifs.csv`: Filtered significant motifs

### Parameters
**Required**:
- `project_name` (str)
- `groupings` (LatchFile): JSON with groupA/groupB cell barcodes
- `archrproject` (LatchDir): Path ending in "_ArchRProject"
- `genome` (Enum): "hg38" or "mm10"

### Creating compare_config.json

- Inspect the user’s input folder and automatically load any `combined_sm_ge.h5ad` file found into an AnnData object.

- Visualize H5AD file with `w_h5`

```python
import json
from pathlib import Path

# Infer comparison column from user request
group_column = "condition"  # or "sample", "cluster"
group_a_value = "Cirrhotic"
group_b_value = "Healthy"

# Get barcodes for each group
group_a_bcs = adata.obs[adata.obs[group_column] == group_a_value].index.tolist()
group_b_bcs = adata.obs[adata.obs[group_column] == group_b_value].index.tolist()

groupings = {
    "groupA": group_a_bcs,
    "groupB": group_b_bcs
}

# Save locally
local_path = "compare_config.json"
with open(local_path, "w") as f:
    json.dump(groupings, f)

# Upload to Latch Data
remote_path = "latch:///compare_config.json"
latch_path = LPath.upload(Path(local_path), remote_path)

# Wrap as LatchFile for workflow input
groupings_file = LatchFile(remote_path)
```

### Formatting Correct Input to `archrproject`

- When using `w_ldata_picker` to populate the `ArchRProject` path, always extract the LData path string via `widget.value.path` before passing it into LatchDir(...).

- Automatically search the user’s input for the folder that ends with `_ArchRProject`. Only ask the user to specify it if you cannot find one.

### Launch Comparison and Wait for it to Complete
```python
params = {
    "project_name": "my_comparison",
    "groupings": groupings_file,
    "archrproject": LatchDir("latch:///Kostallari_SOW313_ATAC_ArchRProject"),
    "genome": "hg38",
}

w = w_workflow(
    wf_name="wf.__init__.compare_workflow",
    key="comparison_workflow_run_1",
    version="0.7.1-8484d6-wip-4ae938",
    params=params,
    automatic=True,
    label="Launch Comparison Workflow"
)

execution = w.value

if execution is not None:
  res = await execution.wait()

  if res is not None and res.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
      # inspect workflow outputs for downstream analysis
      workflow_outputs = list(res.output.values())
```
