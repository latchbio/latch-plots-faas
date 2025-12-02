## Takara Seeker and Trekker Analysis Workflow

Make sure you understand if the kit is Seeker 3x3, Seeker 10x10 or Trekker. If
this information was not provided, ask the user before proceeding.

<plan>
1. **Data Loading** -> `takara/data_loading.md`
2. **Background Removal** (*Seeker only.*) -> `takara/background_removal.md`
3. **Quality Control + Filtering** -> `takara/qc.md`
5. **Normalization** -> `takara/normalization.md`
6. **Feature Selection** -> `takara/feature_selection.md`
7. **Dimensionality Reduction** -> `takara/dimensionality_reduction.md`
8. **Clustering** -> `takara/clustering.md`
9. **Differential Gene Expression (DGE)** -> `takara/diff_gene_expression.md`
10. **Cell Type Annotation** -> `takara/cell_typing.md`
</plan>

<pre_analysis_step>

MANDATORY: Invoke the `redeem_package` tool to install required Takara tools into the workspace.
  - `package_code`: `3015c6c63ecc3f2cd410ea340a36af05777`
  - `package_version_id`: `192`

</pre_analysis_step>

## Launch Trekker Workflow

If user provides FastQ files for the Trekker experiment, launch the Trekker workflow using the `w_workflow` widget.
The workflow requires precise user input because each field maps directly to workflow parameters in the code.
You must **parse user answers**, **normalize them into the required formats**, and then construct the `params` dictionary exactly as shown in the example.

### FASTQ Files

- Every Trekker sample has **two FASTQ files** (paired-end sequencing).
  - **Read 1** → maps to workflow param `fastq_cb`
  - **Read 2** → maps to workflow param `fastq_tags`

These must be provided by the user and wrapped as `LatchFile(latch://...)`.

### Required User Inputs

For each sample, you must ALWAYS **provide a form using latch widgets** for the following values. Each one has a clear mapping to a workflow parameter:

- **Sample ID** → `sample_id`
- **Analysis date (YYYY-MM-DD)** → `analysis_date`
  - Must be normalized to string format `"YYYYMMDD"` when inserted into params.
- **Tile ID** → `tile_id`
- **Single-cell platform** → `sc_platform`
  - User may one of the three options "10x Chromium Next GEM 3’v3.1", “10x Chromium GEM-X 3’v4”, or “BD Rhapsody WTA”; each maps to the following Python strings:
    - “10x Chromium Next GEM 3’v3.1” → `"TrekkerU_C"`
    - “10x Chromium GEM-X 3’v4” → `"TrekkerU_CX"`
    - “BD Rhapsody WTA” → `"TrekkerU_R"`
- **Input directory** (directory on Latch Data that stores the output of the single-cell platform used during the Trekker experiment) → `sc_outdir` (`LatchDir`)
- **Output directory** (where to save Trekker workflow results) → `output_dir` (`LatchDir`)
  - If not provided, default to:
    ```python
    LatchDir("latch://38771.account/Trekker_Outputs/Test_Data")
    ```

When you use the w_ldata_picker widget to populate the `output_dir` or `sc_outdir` params, ALWAYS retrieve the LData path string by accessing the widget `.value.path` before passing to LatchFile(...) or LatchDir(...)

Use the code below as a template, that uses w_workflow. Always use the `automatic` argument or the workflow will not launch. The workflow will launch automatically when the cell is run. Subsequent cell runs with the same key will not relaunch the workflow, so change the key to a new value if you need to relaunch the workflow.

Finally, you need to make sure to wait for the workflow to complete before proceeding. This is included in the code below.

### Example Implementation

```python
from lplots.widgets.workflow import w_workflow
from latch.types import LatchFile, LatchDir

params = {
    "sample_id": "TrekkerU_C_MouseBrain",
    "analysis_date": "20240916",
    "tile_id": "LTTag0053_003",
    "fastq_cb": LatchFile("latch://38771.account/Trekker_Example_Datasets/Test_Data/TrekkerU_C_Mouse_brain_R1_001.fastq.gz"),
    "fastq_tags": LatchFile("latch://38771.account/Trekker_Example_Datasets/Test_Data/TrekkerU_C_Mouse_brain_R2_001.fastq.gz"),
    "sc_outdir": LatchDir("latch://38771.account/Trekker_Example_Datasets/Test_Data/scRNAseq"),
    "sc_platform": "TrekkerU_C",
    "subsample_update": "no",
    "output_dir": LatchDir("latch://38771.account/Trekker_Example_Datasets/Test_Data/")
}

w = w_workflow(
    wf_name="wf.__init__.trekker_pipeline_wf",
    key="trekker_workflow_run_1",
    version="0.2.3-4fdda8",
    params=params,
    automatic=True,
    label="Trekker workflow",
)
execution = w.value

if execution is not None:
  res = await execution.wait()

  if res is not None and res.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
      # inspect workflow outputs for downstream analysis
      workflow_outputs = list(res.output.values())
```
