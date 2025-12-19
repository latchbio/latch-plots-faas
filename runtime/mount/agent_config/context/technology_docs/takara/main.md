## About this doc

- <pre_analysis_questions> any questions to ask *before* analysis if they are not obvious from context 
- <plan> the names of the steps and where to find step docs
- <self_eval_criteria> specific, often numerical, sanity checks after you think you've completed the entire plan

## About the step docs

Each step in the <plan> has its own document you must load before executing the step.

Description of step document tags:

- <goal> describes the scientific goal of the step
- <method> contains a description of the procedure to accomplish goal
- <workflows> contain the names of any Latch workflows you should invoke
- <library> contain the names of any technology specific library could you should use
- <self_eval_criteria> contain specific, often numerical, sanity checks you should run through before determining the step is complete

Make sure you pay close attention to each of these tags when planning, executing and submitting work for each step.

### More information on <workflows>

The value in the tags tells you which workflow document to retrieve in the `wf` dir
nested in the technology dir, eg. `takara/wf`. This document holds info about
parameters, outputs, example usage. What follows is generic information about
how to use workflows:

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
- What tissue and disease conditions describe your data?
- Is the kit type Seeker 3x3, Seeker 10x10 or Trekker?
- Does the H5AD have one or multiple samples?
</pre_analysis_questions>

# TODO: lets put this somewhere else...
<pre_analysis_step>
MANDATORY: Invoke the `redeem_package` tool to install required Takara tools into the workspace.
  - `package_code`: `3015c6c63ecc3f2cd410ea340a36af05777`
  - `package_version_id`: `192`
</pre_analysis_step>

<plan>
1. Reads to Counts (*FastQ ONLY*) -> `steps/reads_to_counts.md`
2. Data Loading -> `steps/data_loading.md`
3. Background Removal (*Seeker ONLY*) -> `steps/background_removal.md`
4. Quality Control + Filtering -> `steps/qc.md`
5. Normalization -> `steps/normalization.md`
6. Feature Selection -> `steps/feature_selection.md`
7. Dimensionality Reduction -> `steps/dimensionality_reduction.md`
8. Clustering -> `steps/clustering.md`
9. Differential Gene Expression  -> `steps/diff_gene_expression.md`
10. Cell Type Annotation -> `steps/cell_typing.md`
</plan>

<self_eval_criteria>
</self_eval_criteria>
