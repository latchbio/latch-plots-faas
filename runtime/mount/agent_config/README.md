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
│   ├── agent_scratch/            # Agent-created notes and state
│   └── notebook_context/         # Auto-generated runtime state
│       ├── cells.md              # (created at runtime)
│       ├── globals.md            # (created at runtime)
│       └── signals.md            # (created at runtime)
└── README.md                     # This file
```

## Modifying Agent Behavior

### 1. Edit System Instructions (`system_prompt.md`)

The `system_prompt.md` file contains all behavior instructions and documentation guidance.

### 2. Add Documentation Files

The agent accesses documentation on-demand using file tools. To add new documentation:

1. **Create the documentation file** in the appropriate subdirectory:
   - Platform workflows → `context/technology_docs/`
   - Latch APIs → `context/latch_api_docs/`

2. **Update `system_prompt.md`** to reference the new documentation:

Add to the `<workflow_intake>` or `<documentation_access>` section:

```markdown
## Technology Workflows

When user mentions a spatial assay platform, read the corresponding workflow:

- **Takara Seeker/Trekker** → `technology_docs/takara_workflow.md`
- **AtlasXOmics** → `technology_docs/atlasxomics.md`
- **Vizgen MERFISH** → `technology_docs/vizgen_workflow.md`
- **Proteomics** → `technology_docs/proteomics_workflow.md`  # ← New addition
```

The agent will use this guidance to decide when to read your documentation.

## How It Works

### Initialization Flow

1. Agent starts up (`agent.py::handle_init`)
2. Agent loads `system_prompt.md` for complete instructions
3. System prompt is loaded directly into the agent
4. Agent is initialized with all behavioral guidance

### On-Demand Documentation Access

When the agent needs specific information (e.g., user mentions "Vizgen"), it:
1. Checks the `<workflow_intake>` or `<documentation_access>` section in its system prompt
2. Uses `read_file` tool to load `technology_docs/vizgen_workflow.md`
3. Follows the workflow steps from that documentation

## Testing Changes

1. Changes to `system_prompt.md` apply on next agent initialization
2. Changes to docs in `context/` are available immediately (agent reads on-demand)

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

2. **Update the system prompt:**
   
   Edit `system_prompt.md` in the `<workflow_intake>` section:
   ```markdown
   ## Workflow Documentation
   
   Once identified, read corresponding workflow from `technology_docs/`:
   - Takara → `takara_workflow.md`
   - AtlasXOmics → `atlasxomics.md`
   - Vizgen MERFISH → `vizgen_workflow.md`
   - Proteomics → `proteomics_workflow.md`  # ← Added
   
   Follow the documented workflow steps precisely.
   ```

3. **Optionally add trigger guidance:**
   
   In the `<workflow_intake>` section, you can add platform specific identification guidance:
   ```markdown
   ## Assay Identification
   
   First, identify the spatial assay platform:
   - Takara Seeker/Trekker
   - Visium
   - Xenium
   - MERFISH (Vizgen)
   - AtlasXOmics
   - Proteomics  # ← Added
   - Other platforms
   ```

That's it! The agent will now read the proteomics docs when relevant.
