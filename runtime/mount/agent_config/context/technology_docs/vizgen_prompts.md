### More information on <workflows>

The value in the tags tells you which workflow document to retrieve in the `wf` dir. This document holds info about parameters, outputs, example usage. What follows is generic information about how to use workflows:

#### Parameter construction

- Provide a form using latch widgets for parameter values.

- **parse user answers**, **normalize them into the required formats**, and then construct the `params` dictionary exactly as shown in the example.

- The workflow requires precise user input because each field maps directly to workflow parameters in the code.

- When you use the w_ldata_picker widget to populate file or directory values, ALWAYS retrieve the LData path string by accessing the widget `.value.path` before passing to LatchFile(...) or LatchDir(...)

#### Launching workflow

Use the code below as a template, that uses w_workflow. Always use the `automatic` argument or the workflow will not launch. The workflow will launch automatically when the cell is run. Subsequent cell runs with the same key will not relaunch the workflow, so change the key to a new value if you need to relaunch the workflow.

Finally, you need to make sure to wait for the workflow to complete before proceeding. This is included in the code below.

### More information on <library>

---

<pre_analysis_questions>

- Is the data already segmented, or do you need to run cell segmentation?

- Does the H5AD have one or multiple samples/regions?

</pre_analysis_questions>

<plan>

1. Cell Segmentation (*if needed*) -> `steps/cell_segmentation.md`

2. Data Loading -> `steps/data_loading.md`

3. Preprocessing -> `steps/preprocessing.md`

4. Quality Control + Filtering -> `steps/qc.md`

5. Spatial Analysis -> `steps/spatial_analysis.md`

6. Secondary Analysis (Cell Type Annotation + Domain Detection) -> `steps/vizgen_secondary_analysis.md`

</plan>

<self_eval_criteria>

</self_eval_criteria>

