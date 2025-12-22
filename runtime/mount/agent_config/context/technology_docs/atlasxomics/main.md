## About this doc

- <pre_analysis_questions> any questions to ask *before* analysis if they are not obvious from context.
- <plan> the names of the steps and where to find step docs
- <data_structure> the organization of AtlasXomics data in the customer's workspace
- <self_eval_criteria> specific, often numerical, sanity checks after you think you've completed the entire plan

## About the step docs

Each step in the <plan> has its own document you must load before executing the step.

Description of step document tags:

- <goal> describes the scientific goal of the step
- <method> contains a description of the procedure to accomplish goal
- <workflows> contain the names of any Latch workflows you should invoke
- <library> contain the names of any technology specific library you should use
- <self_eval_criteria> contain specific, often numerical, sanity checks you should run through before determining the step is complete

Make sure you pay close attention to each of these tags when planning, executing and submitting work for each step.

### More information on <workflows>

The value in the tags tells you which workflow document to retrieve in the `wf` dir
nested in the technology dir, eg. `atlasxomics/wf`. This document holds info about
parameters, outputs, example usage. What follows is generic information about
how to use workflows:

#### Parameter construction

- If a step requires launching a workflow but required inputs are missing, generate a form latch widgets to collect user inputs.
- **parse user answers**, **normalize them into the required formats**, and then construct the `params` dictionary exactly as shown in the example.
- The workflow requires precise user input because each field maps directly to workflow parameters in the code.
- When you use the w_ldata_picker widget to populate file or directory values, ALWAYS retrieve the LData path string by accessing the widget `.value.path` before passing to LatchFile(...) or LatchDir(...)

#### Launching workflow

Use the code below as a template, that uses w_workflow. Always use the `automatic` argument or the workflow will not launch. The workflow will launch automatically when the cell is run. Subsequent cell runs with the same key will not relaunch the workflow, so change the key to a new value if you need to relaunch the workflow.
Finally, you need to make sure to wait for the workflow to complete before proceeding. This is included in the code below.

---

<pre_analysis_questions>
- What organism is your data from (Human - hg38, Mouse - mm10, Rat - rnor6)?
- What tissue and experimental conditions describe your data?
- Do you have raw fragment files and spatial directories, or pre-processed H5AD?
</pre_analysis_questions>

<pre_analysis_step>

MANDATORY: Invoke the `redeem_package` tool to install required AtlasXomics tools into the workspace.
  - `package_code`: `2428814b149447a4c354b3cb4520095b77955bf99cb3eedfef20b920a2a7d3d7`
  - `package_version_id`: `405`

</pre_analysis_step>

<plan>
1. Quality Control + Filtering -> `steps/qc.md`
2. Clustering -> `steps/clustering.md`
3. Differential Analysis -> `steps/de.md`
4. Cell Type Annotation -> `steps/cell_type_annotation/overview.md`
</plan>

<self_eval_criteria>
</self_eval_criteria>

<data_structure>

### Raw Data Paths
**Internal Workspace (13502)**:
- Fragments: `/chromap_outs/[Run_ID]/chromap_output/fragments.tsv.gz`
- Spatial: `/Images_spatial/[Run_ID]/spatial`
- Downstream-analysis ready: `/snap_outs/[project_name]/`

**Collaborator Workspaces**:
- Fragments: `.../Raw_Data/[Run_ID]/chromap_output/fragments.tsv.gz`
- Spatial: `.../Raw_Data/[Run_ID]/spatial`
- Downstream-analysis ready: `.../Processed_Data/[project_name]`

**Key Analysis Files**
- `combined_sm_ge.h5ad`: **Gene activity scores** for all spots/cells across all samples. **Recommended for analyses**
  - .X matrix: gene activity (imputed from chromatin accessibility)
  - Used for: gene expression analysis, cell type annotation
- `combined_sm_motifs.h5ad`: **Motif enrichment scores** for all spots/cells
  - .X matrix: TF motif enrichment scores (870 motifs)
  - Used for: transcription factor activity analysis
- `*_ArchRProject/`: ArchR project directory for R-based analysis

Additional outputs:
- `combined_ge.h5ad`: Gene activity without smoothing
- `combined.h5ad`: Original peak/tile matrix
- `[sample]_g_converted.h5ad`: Per-sample gene activity
- `[sample]_m_converted.h5ad`: Per-sample motif enrichment
- `compare_config.json`: Example grouping file for comparisons
- `cluster_coverages/`, `condition_coverages/`, `sample_coverages/`: BigWig coverage tracks
- `figures/`: QC and analysis plots
- `tables/`: Summary statistics

</data_structure>
