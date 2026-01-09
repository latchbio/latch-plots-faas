# Step 0 — Data Preparation / Data Transform

<goal>
Convert raw Xenium output directories into an `.h5ad` plus viewer-ready assets (e.g. pmtiles) for downstream analysis.
</goal>

<method>
Auto-detect raw Xenium outputs and launch `xenium_preprocess_workflow` via `w_workflow` when no `.h5ad` is attached.
</method>

<workflows>
- `wf.__init__.xenium_preprocess_workflow`

```python
from lplots.widgets.workflow import w_workflow
from latch.types import LatchDir

params = {
    "input_file": LatchDir("latch://38438.account/Scratch/xenium/Input/Xenium_V1_FFPE_TgCRND8_17_9_months_outs"),
    "run_name": "my_run",
    "output_directory": LatchDir("latch:///Xenium_Preprocessing"),
}

w = w_workflow(
    wf_name="wf.__init__.xenium_preprocess_workflow",
    version=None,
    label="Launch Data Preparation Workflow",
    params=params,
    automatic=True,
)
execution = w.value
if execution is not None:
    res = await execution.wait()

    if res is not None and res.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
        workflow_outputs = list(res.output.values())
```
</workflows>

<library>
- `lplots.widgets.workflow`
- `latch.ldata.path`
</library>

<self_eval_criteria>
- All expected Xenium raw files are present in the input directory.
- When this step is run, the attached Xenium output directory contains the following files:
    - `morphology_focus.ome.tif`
    - `morphology_mip.ome.tif`
    - `analysis.tar.gz`
    - `transcripts.parquet`
    - `cell_boundaries.parquet`
    - `cell_feature_matrix.h5`
    - `cells.parquet`
- Workflow completes and produces a processed `.h5ad` and viewer assets.
</self_eval_criteria>
