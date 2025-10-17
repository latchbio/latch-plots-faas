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

    judge_text = """
    Evaluate whether the agent can perform promoter annotation on spatial ATAC-seq peaks.

    Check:
    1. Load spatial ATAC-seq peaks.
    2. Define promoter regions as ±1 kb around the transcription start site (TSS) of genes.
    3. Identify which peaks overlap promoter regions.
    4. Return a list or count of promoter-associated peaks.
    5. Provide at least one example of a peak-gene assignment (e.g. “peak overlaps promoter of ERBB4”).

    """

    create_json(
        "promoter_annotation.json",
        workflow_id="promoter_annotation_001",
        task="Annotate promoter-associated peaks (TSS ±1 kb)",
        data_node="latch://38438.account/Human_DRG_Price_plotsLite/Human_DRG_Price_plotsLite_ArchRProject",
        judge_prompt=judge_text
    )
