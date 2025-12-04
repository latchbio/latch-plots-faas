## Seeker Background Removal Workflow

`seeker_background_removal` removes off-tissue background beads using a three-step density-based filtering algorithm.

*IMPORTANT*: Only use for Curio Seeker data.

### Determining the `min_log10_umi` Threshold

- Plot a histogram of `log10(UMI)` to identify the local minimum between background and tissue signal peaks.

### Determining Kit Type

- Use `'TEN_BY_TEN'` or `'THREE_BY_THREE'` depending on the
  information provided in the data loading step. Ask the user if you still
  don't know and forgot to ask earlier.

### Launch Background Removal and Wait for Completion

```python
from latch.types.file import LatchFile
from latch.types.directory import LatchOutputDir

local_path = "sample.h5ad"
remote_path = "latch:///seeker_data/sample.h5ad"
latch_path = LPath.upload(Path(local_path), remote_path)

params = {
    "input_file": LatchFile("latch:///seeker_data/sample.h5ad"),
    "output_directory": LatchOutputDir("latch:///seeker_data/background_removed/"),
    "sample_id": "my_sample",
    "kit_type": 'TEN_BY_TEN',
    "min_log10_umi": 1.4,  # Adjust based on histogram
    "m": 40,
    "n": 100,
    "p": 5,
    "q": 10,
}

w = w_workflow(
    wf_name="wf.__init__.seeker_background_removal",
    key="background_removal_run_1",
    version="0.1.0-385af1",
    params=params,
    automatic=True,
    label="Launch Background Removal"
)

execution = w.value

if execution is not None:
    res = await execution.wait()

    if res is not None and res.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
        workflow_outputs = list(res.output.values())
```

### Outputs

Look for the `{sample_id}*_background_removed.h5ad` in the output directory.
