# Latch Plots Agent

This directory contains agent context files which can be modified to customize the agent's behavior.

## Overview
- **`system_prompt.md`** - Agent system instructions including documentation guidance
- **`context/`** - Documentation library accessed via file tools

## File Structure

```
agent_config/
├── system_prompt.md              # System prompt
├── context/
│   ├── technology_docs/          # Platform-specific documentation
│   │   ├── vizgen.md
│   │   ├── atlasxomics.md
│   │   └── takara.md
│   ├── latch_api_docs/           # Latch API references
│   │   ├── lpath.md
│   │   └── plots_docs/
│   │       ├── custom-plots.mdx
│   │       ├── reactivity.mdx
│   │       └── widget-types.mdx
│   ├── agent_scratch/            # Agent-created notes and state
│   └── notebook_context/         # Auto-generated runtime state
│       ├── cells.md              # (populated by the agent calling tool to refresh the context)
│       ├── globals.md            # (populated by the agent calling tool to refresh the context)
│       └── signals.md            # (populated by the agent calling tool to refresh the context)
└── README.md                     # This file
```

## Modifying Agent Behavior

### 1. Edit System Instructions (`system_prompt.md`)

The `system_prompt.md` file contains all behavior instructions and documentation guidance.

### 2. Add Documentation Files

The agent accesses documentation on-demand using file tools. To add new documentation:

1. **Create the documentation file** in the appropriate subdirectory:
   - Assay platform documentation → `context/technology_docs/`
   - Latch APIs → `context/latch_api_docs/`

2. **Update `system_prompt.md`** to reference the new documentation:

Add to the `<workflow_intake>` and `<documentation_access>` section:

```markdown
## Assay Platform Documentation

Once identified, read corresponding documentation from `technology_docs/`:
- Takara → `takara.md`
- AtlasXOmics → `atlasxomics.md`
- Vizgen MERFISH → `vizgen.md`
- Proteomics → `proteomics.md`  # ← New addition
```

The agent will use this guidance to decide when to read your documentation.

## How It Works

### Initialization Flow

1. Agent starts up (`agent.py::handle_init`)
2. Agent loads `system_prompt.md` for complete instructions
3. System prompt is loaded directly into the agent

### On-Demand Documentation Access

When the agent needs specific information (e.g., user mentions "Vizgen"), it:
1. Checks the `<workflow_intake>` or `<documentation_access>` section in its system prompt
2. Uses `read_file` tool to load `technology_docs/vizgen.md`
3. Follows the documented steps from that file

## Testing Changes

1. Changes to `system_prompt.md` apply on next agent initialization
2. Changes to docs in `context/` are available immediately (agent reads on-demand)

## Examples

### Example: Add Proteomics Support

1. **Create assay platform documentation:**

   Create `context/technology_docs/proteomics.md`:
   ```markdown
   ## Proteomics Analysis Workflow
   
   1. **Data Loading** - Load mass spec data using Pandas
   2. **Quality Control** - Filter low-quality peptides
   3. **Normalization** - Apply log transformation
   4. **Statistical Analysis** - Run differential expression tests
   ...
   ```

2. **Update the system prompt:**
   
   Edit `system_prompt.md` in the `<documentation_access>` section

3. **Optionally add trigger guidance:**
   
   In the `<workflow_intake>` section, you can add platform specific identification guidance.
