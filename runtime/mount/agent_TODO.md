# Agent Improvements TODO

## Type Safety for Result Handling

Currently using `hasattr()` checks in `handle_query()` for result processing. Should be replaced with proper typing once functionality is confirmed.

### Current code (lines ~453-474):
```python
if AGENT_DEBUG and hasattr(result, "new_items"):
    for item in result.new_items:
        if hasattr(item, "raw_item") and hasattr(item.raw_item, "reasoning"):
            reasoning = item.raw_item.reasoning
            if reasoning and hasattr(reasoning, "content"):
                print(f"[reasoning] {reasoning.content}")

response_content = ""
structured_output = None

if hasattr(result, "final_output_as"):
    try:
        structured_output = result.final_output_as(NotebookResponse)
    except Exception as e:
        print(f"[agent] Could not extract structured output: {e}")

if hasattr(result, "content"):
    response_content = str(result.content)
elif hasattr(result, "output"):
    response_content = str(result.output)
else:
    response_content = str(result)
```

### Should be replaced with:
```python
from agents.result import RunResult
from agents.items import ReasoningItem

# Type the result properly
result: RunResult = await Runner.run(...)

# Use isinstance checks instead of hasattr
if AGENT_DEBUG:
    for item in result.new_items:
        if isinstance(item, ReasoningItem):
            reasoning = item.raw_item.reasoning
            if reasoning and reasoning.content:
                print(f"[reasoning] {reasoning.content}")

# Use typed access
structured_output = result.final_output_as(NotebookResponse)
response_content = str(result.final_output)
```

This will provide:
- Better IDE autocomplete
- Compile-time type checking
- Clearer code intent
- No runtime hasattr overhead