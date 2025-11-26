This document provides ultimate guideline to perform secondary analysis on vizgen MERFISH datasets. Be honest and follow this document to the word.

### ** Cell Type Annotation**
- Prompt the user for **tissue** and **organism**.  
- Assign biological meaning to clusters using **marker-gene dictionaries**.  
- Display **dot plots** of cell types and top markers.

## ** Spatial Domain Detection and Niche Identification

- Rename `adata.obsm["Spatial"]` to `adata.obsm["X_spatial"]`, then save the processed object to `.h5ad`.
    ```python
    if 'Spatial' in adata.obsm.keys():
        adata.obsm['X_spatial'] = adata.obsm['Spatial']

    # Save the processed H5AD file
    output_h5ad_path = "/tmp/merfish_processed_for_domain_detection.h5ad"
    adata.write_h5ad(output_h5ad_path)
    output_lpath = LPath("latch:///merfish/analysis_output/merfish_processed_for_domain_detection.h5ad")
    output_lpath.upload_from(output_h5ad_path)
    local_path = Path(output_h5ad_path)
    output_lpath.upload_from(local_path)
    ```
- Allow the user to a filename and a location to save the input file for the workflow
- Use this file as input for the domain-detection workflow.
- Allow the user to select lambda, runname and output_dir using a form created with **lplots widgets**
```python
params = {
    "input_file": LatchFile("latch:///merfish_processed_for_domain_detection.h5ad"),
    "run_name": "my_run",
    "output_dir": LatchDir("latch:///Domain_detection_output"),
    "lambda_list": "0.5",
}
w = w_workflow(
    wf_name="wf.__init__.domain_detection_wf",
    key="domain_detection_workflow_run_1",
    version=None,
    automatic=True,
    label="Launch Domain Detection Workflow",
    params=params,
)
execution = w.value

if execution is not None:
  res = await execution.wait()

  if res is not None and res.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
      # inspect workflow outputs for downstream analysis
      workflow_outputs = list(res.output.values())
```
