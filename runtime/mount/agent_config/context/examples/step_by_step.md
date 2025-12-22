## Example: Cautious Execution (Step-by-Step Mode)

**User Request**: "Filter this dataset."

**Turn 1: One Step at a Time**
Agent plans but asks for confirmation before acting.
```python
submit_response(
    plan=[{"id": "filter", "description": "Filter cells < 200 genes", "status": "todo"}],
    summary="I plan to filter cells with fewer than 200 genes.",
    questions="Does this threshold work for you?",
    continue=False,
    next_status="done"
)
```

**Turn 2: User Confirms -> Atomic Execution**
User: "Yes, proceed."
Agent creates the filter cell and asks for feedback.

**Example execution cells (step-by-step compliant)**

This demonstrates the constraint-friendly pattern:
- Cell A selects a file via `w_ldata_picker`, downloads via cached `LPath.download`, and **writes** a `Signal`
- Cell B **subscribes** to that signal, visualizes evidence with `w_plot`, and shows filtered results with `w_table`

```python
# Cell A — Select remote file, download locally (cached), write to Signal (no subscription)
from lplots.widgets.ldata import w_ldata_picker
from lplots.widgets.text import w_text_output
from lplots.reactive import Signal
from latch.ldata.path import LPath
from pathlib import Path

if "local_csv_path" not in locals():
    local_csv_path = Signal(None)  # str | None
picked = w_ldata_picker(label="Dataset (CSV) in Latch Data", file_type="file")
lp = picked.value
if lp is None:
    w_text_output(content="Pick a CSV file to continue.", appearance={"message_box": "info"})
else:
    lp: LPath = lp  # already an LPath; use directly
    node_id = lp.node_id()
    if node_id is None:
        w_text_output(content="Could not resolve node_id() for selection.", appearance={"message_box": "danger"})
    else:
        local_p = Path(f"{node_id}{Path(lp.name() or '').suffix}")
        lp.download(local_p, cache=True)
        local_csv_path(local_p.as_posix())  # write-only in this cell
        w_text_output(content=f"Cached download: `{local_p}`", appearance={"message_box": "success"})
```

```python
# Cell B — Subscribe to Signal, show evidence, filter by a threshold, and render via widgets
from lplots.widgets.text import w_text_output, w_text_input
from lplots.widgets.table import w_table
from lplots.widgets.plot import w_plot
from lplots.reactive import Signal
from pathlib import Path
import pandas as pd
import plotly.express as px

if "local_csv_path" not in locals():
    local_csv_path = Signal(None)
local_p_str = local_csv_path()  # subscribes this cell
if local_p_str is None:
    w_text_output(content="No file downloaded yet—run the previous cell first.", appearance={"message_box": "warning"})
else:
    df_raw = pd.read_csv(Path(local_p_str))
    num_cols = df_raw.select_dtypes(include="number").columns
    if len(num_cols) == 0:
        w_text_output(content="No numeric columns found; cannot demonstrate threshold filtering.", appearance={"message_box": "info"})
    else:
        col = num_cols[0]
        thresh_str = w_text_input(label=f"Keep rows where {col} ≥", default="200").value
        try:
            threshold = float(thresh_str)
        except Exception:
            w_text_output(content="Threshold must be numeric (e.g., 200).", appearance={"message_box": "danger"})
        else:
            fig_hist = px.histogram(df_raw, x=col, nbins=50, title=f"Distribution of {col} (choose threshold step-by-step)")
            w_plot(label="Evidence: distribution before filtering", source=fig_hist)
            df_filtered = df_raw[df_raw[col] >= threshold].copy()
            w_text_output(
                content=f"Filtered rows: {df_raw.shape[0]:,} → {df_filtered.shape[0]:,} (threshold={threshold:g}).",
                appearance={"message_box": "primary"},
            )
            # NOTE: `w_table` expects a named global variable (avoid passing expressions like df.head() directly)
            df_filtered_preview_50 = df_filtered.head(50)
            w_table(label="Filtered preview (first 50 rows)", source=df_filtered_preview_50)
```

```python
# [Agent creates the cells above, then pauses for review]
submit_response(
    plan=[{"id": "filter", "description": "Filter cells < 200 genes", "status": "in_progress"}],
    summary="I have created the filter cell.",
    questions="Do the results look fine to you, or would you like to adjust the threshold?",
    continue=False,
    next_status="done"
)
```

**Turn 3: Handling "Bad" Results (No Auto-Correct)**
User: "It filtered too many cells! Only 10 left."
Agent does NOT auto-fix. It asks.
```python
submit_response(
    summary="Observation: Retention is very low (10 cells). The 200 gene threshold may be too strict.",
    questions="How would you like to adjust the threshold? (e.g., try 100 genes?)",
    continue=False,
    next_status="done"
)
```
