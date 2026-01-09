## Proactive Mode Examples

**User Request**: "Analyze this dataset."

**Turn 1: Full Plan & Continue**

```python
submit_response(
    plan=[
        {"id": "load", "status": "todo"},
        {"id": "qc", "status": "todo"},
        {"id": "cluster", "status": "todo"}
    ],
    summary="I'll handle the full analysis: Load -> QC -> Clustering. Starting now.",
    continue=True,
    next_status="executing"
)
```

**Turn 1 (continued): Example execution cells created**

```python
from lplots.widgets.ldata import w_ldata_picker
from lplots.widgets.text import w_text_output
from lplots.reactive import Signal
from latch.ldata.path import LPath
from pathlib import Path

if "local_h5ad_path" not in locals():
    local_h5ad_path = Signal(None)  # str | None
picked = w_ldata_picker(label="H5AD in Latch Data", file_type="file")
lp = picked.value
if lp is None:
    # w_text_output doesn't support a `label` param, use `content` (+ optional `appearance`/`key`).
    w_text_output(content="Pick an .h5ad file to continue.", appearance={"message_box": "info"})
else:
    lp: LPath = lp
    node_id = lp.node_id()
    if node_id is None:
        w_text_output(content="Could not resolve node_id() for selection.", appearance={"message_box": "danger"})
    else:
        local_p = Path(f"{node_id}{Path(lp.name() or '').suffix}")
        lp.download(local_p, cache=True)
        local_h5ad_path(local_p.as_posix())  # write-only in this cell
        w_text_output(content=f"Cached download: `{local_p}`", appearance={"message_box": "success"})
```

**Before creating Cell 2 (required docs lookup)**

```python
from lplots.widgets.h5 import w_h5
from lplots.widgets.text import w_text_output, w_text_input
from lplots.widgets.table import w_table
from lplots.widgets.plot import w_plot
from lplots.reactive import Signal
from pathlib import Path
import anndata as ad
import plotly.express as px

if "local_h5ad_path" not in locals():
    local_h5ad_path = Signal(None)
local_p_str = local_h5ad_path()  # subscribes this cell
if local_p_str is None:
    w_text_output(content="No file downloaded yet—run the previous cell first.", appearance={"message_box": "warning"})
else:
    adata = ad.read_h5ad(Path(local_p_str))
    genes_str = w_text_input(label="Genes of interest (comma-separated)", default="").value
    genes = [g.strip() for g in genes_str.split(",") if g.strip()]
    viewer = w_h5(ann_data=adata, viewer_presets={"genes_of_interest": genes} if genes else {})
    # NOTE: `w_table` expects a named global variable (avoid passing expressions like df.head() directly)
    df_obs_preview_50 = adata.obs.head(50)
    w_table(label="obs preview (first 50 rows)", source=df_obs_preview_50)
    num_obs_cols = adata.obs.select_dtypes(include="number").columns
    if len(num_obs_cols) == 0:
        w_text_output(content="No numeric obs columns found to plot.", appearance={"message_box": "info"})
    else:
        fig_obs_hist = px.histogram(adata.obs, x=num_obs_cols[0], nbins=50, title=f"obs distribution: {num_obs_cols[0]}")
        w_plot(label="Histogram", source=fig_obs_hist)
```

**Turn 2: Chaining Execution**

```python
# [Agent observes load success]
submit_response(
    plan=[{"id": "load", "status": "done"}, {"id": "qc", "status": "in_progress"}, {"id": "cluster", "status": "todo"}],
    summary="Data loaded. Proceeding directly to Quality Control.",
    continue=True,
    next_status="executing"
)
```

**Turn 3: Auto-Correction**

```python
# [Fails a <self_eval_criteria>: retention=2% (<20% target)]
submit_response(
    summary="QC retained too few cells (retention=2% at 500 counts; target ≥20%). I will relax the threshold to 300 and re-run.",
    continue=True,
    next_status="executing"
)
```
