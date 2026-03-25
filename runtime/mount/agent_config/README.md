# Latch Plots Agent

This directory contains agent context files which can be modified to customize the agent's behavior.

## Overview
- **`system_prompt.md`** - Agent system instructions including documentation guidance
- **`context/`** - Documentation library accessed via file tools
- **`.claude/skills/`** - Skills auto-discovered by the agent (technology and Latch integration)

## File Structure

```
agent_config/
├── system_prompt.md              # System prompt
├── context/
│   ├── turn_behavior/            # Per-behavior-mode turn guidelines
│   │   ├── proactive.md
│   │   └── step_by_step.md
│   ├── examples/                 # Per-behavior-mode examples
│   │   ├── proactive.md
│   │   └── step_by_step.md
│   ├── latch_api_docs/           # Latch API references
│   │   ├── latch_api_reference.md
│   │   └── plots_docs/
│   │       ├── custom-plots.mdx
│   │       ├── reactivity.mdx
│   │       └── widget-types.mdx
│   └── agent_scratch/            # Agent-created notes and state
└── README.md                     # This file

.claude/skills/                   # Cloned at runtime into /opt/latch/plots-faas/.claude/skills/
├── takara-devkit/                # Takara Seeker/Trekker technology skill
├── xenium-devkit/                # 10X Xenium technology skill
├── vizgen-devkit/                # Vizgen MERFISH technology skill
├── atlasx-devkit/                # AtlasXomics technology skill
├── latch-workflows/              # Latch workflow launching skill
├── latch-plots-ui/               # Latch UI and plot widget skill
└── latch-data-access/            # Latch Data path handling skill
```

## Modifying Agent Behavior

### 1. Edit System Instructions (`system_prompt.md`)

The `system_prompt.md` file contains all behavior instructions and documentation guidance.

### 2. Skills

Skills are auto-discovered from `.claude/skills/` at runtime. Each skill repo must have a `SKILL.md` at the root with YAML frontmatter (`name` and `description`).

There are two kinds of skills:

- **Technology skills** hold platform-specific scientific guidance (workflow order, step docs, helper libraries). They are portable and do not assume a Latch-specific system prompt.
- **`latch-*` skills** hold reusable Latch-specific workflow, UI, and data-access guidance.

Technology skills should prefer matching `latch-*` skills when available but fall back to their own local docs when not.

### 3. Add Documentation Files

The agent accesses documentation on-demand using file tools. To add new documentation:

1. **Create the documentation file** in the appropriate subdirectory under `context/`
2. **Update `system_prompt.md`** to reference the new documentation

## How It Works

### Initialization Flow

1. Agent starts up (`agent.py::handle_init`)
2. Agent loads `system_prompt.md` for complete instructions
3. Skills are auto-loaded from `.claude/skills/` based on user request matching
4. System prompt is loaded directly into the agent

### Skill Loading

When the user mentions a supported platform (e.g., "Takara", "Visium"), the agent:
1. Matches the request against skill descriptions in `.claude/skills/*/SKILL.md`
2. Loads the matching technology skill
3. Follows the skill's workflow order, step docs, and helper-library guidance
4. Uses `latch-*` skills for Latch-specific execution details when available
