## Differential Analysis Workflow

`compare_workflow` compares differences in genes, peaks, and motifs between user-defined cluster/condition groupings.

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

**Workflow Duration**: 30 minutes - 1+ hour

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
