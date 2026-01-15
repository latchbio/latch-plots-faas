## Step-by-Step Mode Examples

**User Request**: "Analyze this dataset."

**Turn 1: Plan Confirmation**

```python
submit_response(
    plan=[
        {"id": "filter", "description": "Filter cells by gene count", "status": "todo"},
        {"id": "qc", "description": "Quality control", "status": "todo"},
        # ...
    ],
    summary="I plan to filter cells with fewer than 200 genes",
    questions="Does this threshold work for you?",
    continue=False,
    next_status="done"
)
```

**Turn 2: User Confirms -> Step Execution**

User: "Yes, proceed."

**Example cells**

```python
from lplots.widgets.ldata import w_ldata_picker
from lplots.widgets.text import w_text_output
from lplots.reactive import Signal
from latch.ldata.path import LPath
from pathlib import Path

if "local_h5ad_path" not in locals():
    local_h5ad_path = Signal(None)  # str | None
picked = w_ldata_picker(label="Dataset (.h5ad) in Latch Data", file_type="file")
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
    num_obs_cols = adata.obs.select_dtypes(include="number").columns
    if len(num_obs_cols) == 0:
        w_text_output(content="No numeric obs columns found; cannot demonstrate threshold filtering.", appearance={"message_box": "info"})
    else:
        col = num_obs_cols[0]
        thresh_str = w_text_input(label=f"Keep cells where obs[{col}] ≥", default="200").value
        try:
            threshold = float(thresh_str)
        except Exception:
            w_text_output(content="Threshold must be numeric (e.g., 200).", appearance={"message_box": "danger"})
        else:
            fig_hist = px.histogram(adata.obs, x=col, nbins=50, title=f"obs distribution: {col} (choose threshold step-by-step)")
            w_plot(label="Evidence: distribution before filtering", source=fig_hist)
            adata_filtered = adata[adata.obs[col] >= threshold].copy()
            w_text_output(
                content=f"Filtered cells: {adata.n_obs:,} → {adata_filtered.n_obs:,} (threshold={threshold:g}).",
                appearance={"message_box": "primary"},
            )

            viewer = w_h5(ann_data=adata_filtered)

            # NOTE: `w_table` expects a named global variable (avoid passing expressions like df.head() directly)
            df_obs_filtered_preview_50 = adata_filtered.obs.head(50)
            w_table(label="Filtered obs preview (first 50 rows)", source=df_obs_filtered_preview_50)
```

```python
# [Agent runs the cells above, then performs self-eval criteria]
# [Self-eval criterion failed: retention too low (10 cells)]
submit_response(
    plan=[
        {"id": "filter", "description": "Filter cells by gene count", "status": "in_progress"},
        {"id": "qc", "description": "Quality control", "status": "todo"},
        # ...
    ],
    summary="The filtering completed. Retention is very low (10 cells). The 200 gene threshold may be too strict.",
    questions="How would you like to adjust the threshold before rerunning? (e.g., 100 or 150 genes)",
    continue=False,
    next_status="awaiting_user_response"
)
```

**Turn 3: User Adjusts Threshold. Agent Executes Cells Successfully and Asks for User Confirmation**

**User Request**: "Try 100 genes"

```python
# Agent changes threshold to 100 genes, reruns cells, and self-eval criteria passes
submit_response(
    plan=[
        {"id": "filter", "description": "Filter cells by gene count", "status": "in_progress"},
        {"id": "qc", "description": "Quality control", "status": "todo"},
        # ...
    ],
    summary="Filtering step complete. Results: 3,873 / 4,992 spots retained (77.6%). Retention is now acceptable.",
    questions="Are you satisfied with these results? Should I proceed to the normalization step?",
    continue=False,
    next_status="awaiting_user_response"
)
```
