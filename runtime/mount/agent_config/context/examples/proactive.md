## Example: Autonomous Analysis (Proactive Mode)

**User Request**: "Analyze this dataset."

**Turn 1: Full Plan & Immediate Action**
Agent plans the entire analysis and starts immediately.
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

**Turn 1 (continued): Example execution cells created immediately**

Below is an example of the kind of *self-contained* notebook code the agent would create right away to satisfy the core constraints:
- file selection via `w_ldata_picker`
- remote file operations via `LPath` (cached download)
- cross-cell dependency via `Signal` (Cell A writes, Cell B reads)
- user-facing output via widgets (`w_text_output`, `w_table`, `w_plot`)

```python
# Cell A — Select remote file, download locally (cached), write to Signal (no subscription)
from lplots.widgets.ldata import w_ldata_picker
from lplots.widgets.text import w_text_output
from lplots.reactive import Signal
from latch.ldata.path import LPath
from pathlib import Path

if "local_csv_path" not in locals():
    local_csv_path = Signal(None)  # str | None
picked = w_ldata_picker(label="CSV in Latch Data", file_type="file")
lp = picked.value
if lp is None:
    w_text_output(content="Pick a CSV file to continue.", appearance={"message_box": "info"})
else:
    # already an LPath; use directly
    lp: LPath = lp
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
# Cell B — Subscribe to Signal, load data, render via widgets
from lplots.widgets.text import w_text_output
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
    # NOTE: `w_table` expects a named global variable (avoid passing expressions like df.head() directly)
    df_preview_50 = df_raw.head(50)
    w_table(label="Preview (first 50 rows)", source=df_preview_50)
    num_cols = df_raw.select_dtypes(include="number").columns
    if len(num_cols) == 0:
        w_text_output(content="No numeric columns found to plot.", appearance={"message_box": "info"})
    else:
        fig_hist = px.histogram(df_raw, x=num_cols[0], nbins=50, title=f"Distribution of {num_cols[0]}")
        w_plot(label="Histogram", source=fig_hist)
```

**Turn 2: Chaining Execution**
Load cell finishes. Agent immediately starts QC without asking.
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
QC removes too many cells. Agent attempts fix automatically.
```python
# [Agent observes low retention]
submit_response(
    summary="QC retained too few cells at 500 counts. I will automatically relax the threshold to 300 and re-run.",
    continue=True,
    next_status="executing"
)
```
