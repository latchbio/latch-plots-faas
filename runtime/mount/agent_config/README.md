# Latch Plots Agent Dev Kit

This directory contains the agent configuration that can be modified to customize the agent's behavior without touching the core codebase.

## Overview

The agent uses two main configuration files:
- **`prompts.py`** - System instructions and documentation references
- **`docs/`** - Workflow documentation files

## File Structure

```
agent_config/
├── prompts.py              # Agent instructions and doc references
├── docs/
│   └── takara_workflow.md  # Example workflow documentation
└── README.md               # This file
```

## Modifying Agent Behavior

### 1. Edit System Instructions (`prompts.py`)

The `system_instruction` variable contains the core agent prompt. You can modify:

**Planning behavior:**
```python
system_instruction = """
...
**Planning**
- Create a plan only for complex tasks (more than 5 steps)  # ← Modified
- Once plan is approved, proceed without re-confirming
...
"""
```

**Output format:**
```python
system_instruction = """
...
**Output JSON Object**
- Use empty array (`[]`) for no items
- Set `questions: null` unless needing input
...
"""
```

### 2. Add Documentation (`external_docs`)

Add new workflow documentation by adding entries to `external_docs`:

```python
external_docs = [
    # ... existing docs ...
    {
        "name": "proteomics_docs",
        "path": "/opt/latch/latch-plots-faas/runtime/mount/agent_config/docs/proteomics_workflow.md",
        "type": "file",
    },
]
```

Then create the file at `docs/proteomics_workflow.md`:

```markdown
## Proteomics Analysis Workflow

1. **Data Loading** - Load mass spec data
2. **Quality Control** - Filter low-quality peptides
3. **Normalization** - Apply log transformation
...
```

The agent will automatically include this as `<proteomics_docs>...</proteomics_docs>` in its instructions.

## How It Works

### Initialization Flow

1. Agent starts up (`agent.py::handle_init`)
2. `config_loader.py` loads `prompts.py`
3. `build_full_instruction()` assembles the complete prompt:
   - Loads external docs (plots_docs, random_pointers, lpath_docs, takara_docs)
   - Wraps them in XML tags
   - Appends `system_instruction`
   - Injects current notebook context
4. Agent is initialized with combined instructions

### Assembled Instruction Format

```
<plots_docs>
<plots_docs_custom-plots>
...
</plots_docs_custom-plots>
</plots_docs>

<random_pointers>
...
</random_pointers>

<lpath_docs>
...
</lpath_docs>

<takara_docs>
...
</takara_docs>

You are a spatial data analysis agent...
[system_instruction content]
...

<notebook_context>
[current notebook state]
</notebook_context>
```

## Testing Changes

After modifying configuration files, restart the agent process to apply changes:

1. Changes apply on next agent initialization
2. No pod rebuild required
3. Check logs for loading errors: `[config_loader] Warning: ...`

## External Documentation

These paths are managed separately and should not be modified here:

- `/opt/latch/nucleus-llm-inference/prompt_components/plots_docs/` - Plot widget documentation
- `/opt/latch/nucleus-llm-inference/prompt_components/random_pointers/` - General instructions
- `/opt/latch/nucleus-llm-inference/prompt_components/lpath.py` - LPath API docs

## Examples

### Example: Add Proteomics Support

1. Create `docs/proteomics_workflow.md` with workflow steps
2. Add to `external_docs` in `prompts.py`:
   ```python
   {
       "name": "proteomics_docs",
       "path": "/opt/latch/plots-faas/runtime/mount/agent_config/docs/proteomics_workflow.md",
       "type": "file",
   }
   ```
3. Modify `system_instruction` to reference it:
   ```python
   - **For Planning**:
       - Ask user about assay type (spatial, proteomics, etc.)
       - If proteomics, follow <proteomics_docs> workflow
   ```

## Troubleshooting

**"Config module not found"**
- Check file path in error message
- Ensure files are in `/opt/latch/plots-faas/runtime/mount/agent_config/`

**Changes not applying**
- Restart agent process
- Verify no syntax errors in modified files
