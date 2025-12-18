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
```python
# [Agent creates cell]
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
