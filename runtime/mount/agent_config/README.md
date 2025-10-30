# Latch Plots Agent Dev Kit

This directory contains the agent configuration that can be modified to customize the agent's behavior without touching the core codebase.

## Overview

The agent uses on-demand documentation access via file manipulation tools:
- **`prompts.py`** - System instructions
- **`context/`** - Documentation library accessed via file tools
- **`context/index.md`** - Documentation catalog (loaded into system prompt)

## File Structure

```
agent_config/
├── prompts.py                    # Agent system instruction
├── context/
│   ├── index.md                  # Doc catalog (loaded on init)
│   ├── technology_docs/          # Platform-specific workflows
│   │   ├── vizgen_workflow.md
│   │   ├── atlasxomics.md
│   │   └── takara_workflow.md
│   ├── latch_api_docs/           # Latch API references
│   │   ├── lpath.md
│   │   └── plots_docs/
│   │       ├── custom-plots.mdx
│   │       ├── reactivity.mdx
│   │       └── widget-types.mdx
│   └── current_context/          # For agent to maintain state
└── README.md                     # This file
```

## Modifying Agent Behavior

### 1. Edit System Instructions (`prompts.py`)

The `system_instruction` variable contains the core agent prompt. You can modify any section:

**Example - Change planning behavior:**
```python
system_instruction = """
...
**Planning Protocol**
- Create a plan only for complex tasks (more than 5 steps)  # ← Modified
...
"""
```

### 2. Add Documentation Files

The agent accesses documentation on-demand using file tools. To add new documentation:

1. **Create the documentation file** in the appropriate subdirectory:
   - Platform workflows → `context/technology_docs/`
   - Latch APIs → `context/latch_api_docs/`
   - Temporary notes → `context/current_context/`

2. **Update `context/index.md`** to list the new documentation:

```markdown
### New Platform
**File**: `technology_docs/proteomics_workflow.md`
**Read when**: User mentions "proteomics" or mass spectrometry
**Contents**: Step-by-step proteomics analysis pipeline...
```

The agent will use the INDEX to decide when to read your documentation.

## How It Works

### Initialization Flow

1. Agent starts up (`agent.py::handle_init`)
2. Agent loads `prompts.py` for `system_instruction`
3. Agent reads `context/index.md` (documentation catalog)
4. System prompt is assembled:
   ```
   <documentation_index>
   [index.md contents]
   </documentation_index>
   
   [system_instruction]
   ```
5. Agent is initialized with this compact prompt

### On-Demand Documentation Access

The agent has 5 file tools available:

- **`glob_file_search`** - Find files by pattern
- **`grep`** - Search text in files with line numbers  
- **`read_file`** - Read file contents (with offset/limit for large files)
- **`search_replace`** - Edit files (for maintaining state)
- **`bash`** - Execute bash commands

When the agent needs specific information (e.g., user mentions "Vizgen"), it:
1. Checks the INDEX to find relevant documentation
2. Uses `read_file` to load `technology_docs/vizgen_workflow.md`
3. Follows the workflow steps from that documentation

This approach keeps the system prompt small while providing access to extensive documentation.

## Testing Changes

After modifying configuration files, restart the agent process to apply changes:

1. Changes to `prompts.py` or `index.md` apply on next agent initialization
2. Changes to other docs in `context/` are available immediately (agent reads on-demand)
3. No pod rebuild required
4. Check logs for loading errors: `[agent] Warning: ...`

## Examples

### Example: Add Proteomics Support

1. **Create workflow documentation:**
   
   Create `context/technology_docs/proteomics_workflow.md`:
   ```markdown
   ## Proteomics Analysis Workflow
   
   1. **Data Loading** - Load mass spec data using Pandas
   2. **Quality Control** - Filter low-quality peptides
   3. **Normalization** - Apply log transformation
   4. **Statistical Analysis** - Run differential expression tests
   ...
   ```

2. **Update the INDEX:**
   
   Add to `context/index.md`:
   ```markdown
   ### Proteomics
   **File**: `technology_docs/proteomics_workflow.md`
   **Read when**: User mentions "proteomics", "mass spec", or "peptides"
   **Contents**: Complete proteomics analysis pipeline from data loading through differential expression.
   ```

3. **Optionally update system instruction:**
   
   In `prompts.py`, you can add a hint in the **Assay Intake** section:
   ```python
   * First, identify the assay type (spatial transcriptomics, proteomics, etc.)
   * Once identified, read the corresponding workflow from technology_docs/
   ```

That's it! The agent will now read the proteomics docs when relevant.

## Troubleshooting

**Agent not reading documentation**
- Check that index.md properly describes when to read each doc
- Verify file paths in index.md match actual file locations
- Check agent logs for file reading errors

**Changes not applying**
- Restart agent process for `prompts.py` or `index.md` changes
- Documentation files in `context/` are read on-demand (no restart needed)
- Verify no syntax errors in modified files
