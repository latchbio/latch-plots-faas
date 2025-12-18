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
