`wf.__init__.rapids_single-cell_preprocessing` is a workflow to perform single-cell/spatial preprocessing operations using the GPU-accelerated scRapids package.

## Save Current AnnData

Save the current `adata` to `ldata`. Use lplots widgets to let the user select both the save location and file name. After saving, confirm to the user that the file was successfully saved.

**ABSOLUTE RULE**: The workflow contains **rigorously tested, gold-standard SnapATAC2 code internally and should ALWAYS be used OVER writing code from scratch.**

---

### Input Data & Decision Logic

1. **adata** → Use the file you saved to `ldata` in **Save Current Data** (`LatchFile`) and map it to `input_file`.  
2. **Output Directory** → Let the user select where workflow results should be written. Map this to `output_directory` (`LatchOutputDir`).  
3. **Run Name** → Let the user specify the run name. Map this to `run_name`.  
4. **Batch Key** → Always check, via lplots widgets, whether the user wants batch correction.  
   - If **Yes**, allow them to choose batch keys from `adata.obs`.  
   - Only show **categorical** columns.  
   Map this to `batch_key`.  
5. **Preprocessing Check** → If `adata` is already preprocessed, set `skip_qc = True`. Map this to `skip_qc`.  
6. **Clustering Resolutions** → Map to `clustering_resolution` (`typing.List[float]`). Default: `[0.3, 0.5, 0.7, 1.0]`.

---

⚠️ **ALWAYS** make sure the **paths** that you pass to the workflow are valid. Once inputs are gathered, construct a parameter dictionary like:

```python
params = {
    "output_directory": LatchOutputDir(
        "latch:///Rapids_SingleCell_Output"
    ),
    "run_name": "run_1",
    "clustering_resolution": [0.3, 0.5, 0.7, 1.0],
}
```


Then launch the workflow using the **exact code** below in the same cell:

```python
w = w_workflow(
    wf_name="wf.__init__.rapids_single-cell_preprocessing",
    version=None,
    params=params,
    automatic=True,
    label="Run scRapids",
)
execution = w.value

if execution is not None:
    res = await execution.wait()
    if res.status == 'SUCCEEDED':
        w_text_output(
            content="✓ RAPIDS workflow completed successfully!",
            appearance={"message_box": "success"}
        )
```

The wait step is essential—it ensures the workflow completes before proceeding. After completion, always load the new adata from {output_directory}/{run_name} on ldata and never skip this.

Wait for the workflow to finish patiently before moving on! 
