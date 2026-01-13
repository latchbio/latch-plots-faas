<goal>
Determine preprocessing path based on dataset size and route to the appropriate workflow or Scanpy-based preprocessing.
</goal>

<method>
1/ **Check dataset size**  
   - Inspect the AnnData object and determine whether it contains **more than 100,000 cells**.

2/ **If `n_cells > 100,000`**  
   - **Always follow the RAPIDS preprocessing instructions** described in `wf/rapids_wf.md`.  
   - Launch and complete the GPU-accelerated preprocessing workflow before moving forward.  
   - Proceed directly to the **Spatial Analysis Step** afterward.

3/ **If `n_cells ≤ 100,000`**  
   - **Skip RAPIDS**.  
   - **Always follow the Scanpy preprocessing instructions** described in `scanpy_preprocessing.md`.  
   - After completing these steps, proceed to the **Spatial Analysis Step**.
</method>


<workflows>
wf/rapids_wf.md

```python
params = {
    "input_file": saved_h5ad,
    "output_directory": LatchOutputDir("latch:///Rapids_Output"),
    "run_name": "rapids_1",
    "clustering_resolution": [0.3, 0.5, 0.7, 1.0],
}

w = w_workflow(
    wf_name="wf.__init__.rapids_single-cell_preprocessing",
    version=None,
    params=params,
    automatic=True,
)

execution = w.value
if execution:
    res = await execution.wait()
    if res.status == "SUCCEEDED":
        w_text_output(
            content="✓ RAPIDS workflow completed!",
            appearance={"message_box": "success"}
        )
```
</workflows>

<library>
</library>

<self_eval_criteria>
</self_eval_criteria>
