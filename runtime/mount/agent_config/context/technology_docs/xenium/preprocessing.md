# Step 2 — Preprocessing (run only if the user confirms)

<goal>
Create a clean, integrated embedding and cluster structure ready for differential expression and cell-type annotation.
</goal>

<method>

- If `adata.n_obs > 100_000`, you **must** use RAPIDS for preprocessing (see `technology_docs/rapids.md`). Otherwise use Scanpy (see `technology_docs/xenium/scanpy.md`).
- Run QC, normalization, PCA, Harmony (if needed), neighbors, UMAP, Leiden and DGE on the Xenium `AnnData`.
</method>

<workflows>
technology_docs/rapids.md

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
- `scanpy`
- `harmonypy` / `sce.pp.harmony_integrate`
- `rapids_singlecell`
- `plotly`
- `matplotlib`
</library>

<self_eval_criteria>
* QC filters remove low-quality cells without dropping all cells.
* PCA/UMAP embeddings and Leiden clusters are present and usable for downstream steps.
</self_eval_criteria>
