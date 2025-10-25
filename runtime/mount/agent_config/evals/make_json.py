import json
from pathlib import Path

def create_json(
    file_path: str,
    workflow_id: str,
    task: str,
    data_node: str,
    judge_prompt: str,
):
    """
    Creates a JSON file with the given fields.
    judge_prompt can be a normal multiline string.
    """
    data = {
        "id": workflow_id,
        "task": task,
        "data_node": data_node,
        "judge_prompt": judge_prompt.strip()
    }

    # json.dumps will automatically convert newlines to \n
    Path(file_path).write_text(json.dumps(data, indent=4))

if __name__ == "__main__":

    base_dir = "/Users/hannahle/Documents/GitHub/latch-plots-faas/runtime/mount/agent_config/evals/atlasxomics/"

    ### Cell type annotation

    judge_text = """
    Evaluate whether the agent successfully finds and labels hepatic stellate cells. 

    Check:
    1. Load and visualize gene activity score file with `w_h5` widget. 
    2. Ask a question to confirm organism and tissue type.
    2. Create a form with a button to collect cell types and marker genes from users. The form is populated with sensible defaults. 
    3. Perform gene set scoring using `scanpy.tl.score_genes`.
    4. Add new label for cell type to `adata.obs`. 

    """

    create_json(
        base_dir + "hsc_typing.json",
        workflow_id="hsc_typing_1",
        task="Identify hepatic stellate cells",
        data_node="latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC",
        judge_prompt=judge_text
    )
    
    ### Clustering

    judge_text_2 = """
    Evaluate whether the agent launches the correct workflow for clustering.

    Check:
    1. The `genome` parameter must be one of `"hg38"`, `"mm10"`, or `"rnor6"`.
    2. The length of the `runs` list must match the number of samples in the `adata` object (typically under `adata.obs["sample"]`), where sample IDs follow the format `DXXXXX_NGYYYYY`.
    3. The `resolution` value should reflect the desired clustering granularity: higher values yield more clusters, lower values fewer.
    4. All parameters above should be exposed in the workflow form.
    5. If the user requests clustering on a subset of the dataset, the agent must create a subsetted `adata` object, upload it to Latch using the correct APIs.
    6. The adata subset created must contain >0 cells. 
    7. Must use the exact workflow name and version in the `w_workflow` widget and include `w.value` in the code block to render the workflow launch button.
    8. Once workflow finishes, retrieve workflow outputs and make visualization. 
    ```
    w = w_workflow(
        wf_name="wf.__init__.opt_workflow",
        version="0.3.5-9e16e4",
            ....
    )
    w.value
    ```
    """


    create_json(
        base_dir + "clustering.json",
        workflow_id="clustering_1",
        task="Help me find more clusters on the dataset",
        data_node="latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC",
        judge_prompt=judge_text_2
    )

    ### Differential Expression

    judge_text_3 = """
    Evaluate whether the agent launches the comparison workflow to compare between motifs / chromatin accessibility between groups of interest. 
    
    Check: 
    1. Create comparison configs between two comparison groups and upload to Latch Data.
    2. Each group must contain more than 0 cells. 
    3. Display a workflow form and populate it with parameters. 
    4. The ArchR project parameter must point to a file that exists on Latch Data. 
    5. Must use the exact workflow name and version in the `w_workflow` widget and include `w.value` in the code block to render the workflow launch button.
    6. Once workflow finishes, retrieve workflow outputs and make visualization. 
    ```
    w = w_workflow(
    wf_name="wf.__init__.compare_workflow",
    version="0.7.1-8484d6-wip-4ae938",
    ...
    )
    w.value
    ```
    """

    create_json(
        base_dir + "de.json",
        workflow_id="de_1",
        task="Perform differential gene activity comparison between two conditions for cluster 5.",
        data_node="latch://38438.account/AtlasxOmics/Kosta/Kostallari_SOW313_ATAC",
        judge_prompt=judge_text_3
    )

    ### General AnnData operations
    judge_text_4 = """
    Evaluate whether the agent merges `adata.obs` column correctly. 

    Check:
    1. Use a `adata.obs` column that exist. 
    2. Use values that exist within that `adata.obs` column.
    3. Create a new `adata.obs` column for merged values. 

    """

    create_json(
        base_dir + "merge_obs.json",
        workflow_id="merge_obs_1",
        task="Merge cluster 1 and 4",
        data_node="latch://38438.account/AtlasxOmics/Stomach_Case_Study_Plots/combined_sm_ge.h5ad",
        judge_prompt=judge_text_4
    )
