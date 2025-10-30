# Documentation Index

This directory contains documentation that you can read on-demand using file manipulation tools (glob_file_search, grep, read_file). Only read documentation when you actually need it - don't load everything at once.

## Available File Tools

- **glob_file_search**: Find files by pattern (e.g., `*.md`, `vizgen*`)
- **grep**: Search for text patterns in files with line numbers
- **read_file**: Read file contents with optional offset/limit for large files
- **search_replace**: Edit files (useful for maintaining notes in current_context/)
- **bash**: Execute bash commands for exploring the directory

All paths below are relative to `agent_config/context/`.

---

## Technology Workflow Documentation

**Location**: `technology_docs/`

**When to read**: User mentions a specific spatial transcriptomics or spatial genomics platform, or asks for platform-specific analysis workflows.

### Vizgen MERFISH
**File**: `technology_docs/vizgen_workflow.md`
**Read when**: User mentions "Vizgen", "MERFISH", "Merscope", or provides MERFISH spatial data
**Contents**: Authoritative step-by-step pipeline for Vizgen MERFISH experiments including data loading, QC, preprocessing, clustering, differential expression, and cell type annotation. Includes specific guidance on spatial image handling and pmtiles.

### AtlasXOmics
**File**: `technology_docs/atlasxomics.md`
**Read when**: User mentions "AtlasXOmics", "AtlasX", "spatial ATAC-seq", or provides AtlasXOmics data (gene activity scores or motif enrichment)
**Contents**: Step-by-step pipeline for AtlasXOmics spatial ATAC-seq analysis. Includes workflow launch instructions for clustering and differential analysis using w_workflow, gene activity analysis, and cell type annotation.

### Takara Seeker/Trekker
**File**: `technology_docs/takara_workflow.md`
**Read when**: User mentions "Takara", "Seeker", "Trekker", or provides Takara platform data
**Contents**: Authoritative pipeline for Takara Seeker and Trekker experiments with Seeker-specific preprocessing steps. Covers background removal (Seeker only), QC, normalization, feature selection, dimensionality reduction, clustering, and cell type annotation.

---

## Latch API & Platform Documentation

**Location**: `latch_api_docs/`

**When to read**: Need to work with Latch-specific APIs, remote file operations, widgets, or custom plotting.

### LPath API
**File**: `latch_api_docs/lpath.md`
**Read when**: 
- Working with remote `latch://` paths
- Downloading files from Latch Data
- Uploading files to Latch Data
- Need to understand LPath vs pathlib.Path usage
- Errors related to remote file operations

**Contents**: Complete LPath API reference for remote file operations. Includes path construction, file downloads, uploads, directory operations, and common patterns. Critical for any interaction with Latch Data storage.

### Widgets Reference
**File**: `latch_api_docs/plots_docs/widget-types.mdx`
**Read when**: 
- Need to create user input widgets (selects, text inputs, checkboxes, etc.)
- Working with ldata pickers for file selection
- Need to understand widget `.value` and reactivity
- Creating interactive controls for parameters

**Contents**: Comprehensive lplots widget API including all widget types (select, ldata_picker, text_input, checkbox, radio, slider, etc.), their usage patterns, and how to access widget values.

### Custom Plots
**File**: `latch_api_docs/plots_docs/custom-plots.mdx`
**Read when**: 
- Creating matplotlib/seaborn visualizations
- Creating plotly interactive plots
- Need to render figures in the UI
- Understanding `fig*` variable naming convention
- Using `w_plot` widget

**Contents**: Patterns for creating custom plots with matplotlib/seaborn (static) or plotly (interactive). Explains variable naming requirements and how to render plots with w_plot widget.

### Reactivity System
**File**: `latch_api_docs/plots_docs/reactivity.mdx`
**Read when**: 
- Working with Signals for cross-cell dependencies
- Need cells to react to widget value changes
- Creating reactive workflows
- Debugging infinite loops or unexpected re-execution
- Understanding `.value` vs `._signal.sample()`

**Contents**: Latch Plots reactivity system documentation. Covers Signal API, reactive execution, subscription patterns, and how to prevent infinite loops. Critical for multi-cell workflows where cells depend on each other.

---

## Current Context (Runtime State)

**Location**: `current_context/`

**When to use**: These files are automatically refreshed before each turn with the latest notebook state.

### Automatically Generated Files

**`cells.md`** - Notebook cell structure
**Read when**: Need to check what cells exist, find specific code, see cell status, or locate widgets
**Contents**: All notebook cells with searchable markers (CELL_ID, CELL_INDEX, TYPE, STATUS). Code is between CODE_START/CODE_END markers. Optimized for grep searching.
**Refresh**: Automatic before each turn

**`globals.md`** - Global variables
**Read when**: Need to check what variables exist, their types, shapes, or properties  
**Contents**: All global variables with TYPE, SHAPE, COLUMNS, DTYPES, etc. Each variable has its own section starting with `## Variable: name`. Optimized for grep searching.
**Refresh**: Automatic before each turn

**`signals.md`** - Reactivity and dependencies
**Read when**: Need to understand which cells depend on each other or debug reactive execution
**Contents**: Signal dependencies between cells, widget signals, and global variable signals. Shows which cells subscribe to which signals.
**Refresh**: Automatic before each turn

### Additional Uses

You can also create your own files here via search_replace to:
- Maintain a running log of analysis steps
- Keep notes on user preferences or decisions
- Store temporary state across turns

---

## Usage Guidelines

### Strategy for Efficient Documentation Access

1. **Start with this INDEX** - Use it to decide which docs to read
2. **Read selectively** - Only read documentation when you encounter a specific need
3. **Use grep for targeted searches** - If you know what you're looking for, grep is faster than reading entire files
4. **Use offset/limit for large files** - Read files in chunks if they're very long
5. **Cache in memory** - Once you've read a doc, remember its key points for the current session

### Example Workflows

**User asks about Vizgen analysis:**
```
1. Check index.md (already in context)
2. read_file: technology_docs/vizgen_workflow.md
3. Follow the workflow steps
```

**User mentions latch:// path issues:**
```
1. Check index.md - see it's LPath related
2. read_file: latch_api_docs/lpath.md (or use offset/limit if just need specific section)
3. Apply LPath patterns
```

**Need to create a widget but forgot the API:**
```
1. grep: pattern="w_select" path="latch_api_docs/plots_docs/widget-types.mdx"
2. Read the specific section returned by grep
```

**User mentions multiple platforms:**
```
1. Clarify which platform they're using
2. Read only the relevant technology_docs file
```

---

## Quick Reference

| Topic | File | Key Triggers |
|-------|------|--------------|
| **Runtime Context** | | **Auto-refreshed each turn** |
| Notebook cells | current_context/cells.md | Check cells, find code, cell status |
| Global variables | current_context/globals.md | Check variables, types, shapes |
| Reactivity | current_context/signals.md | Understand dependencies |
| **Technology Workflows** | | |
| Vizgen MERFISH | technology_docs/vizgen_workflow.md | "Vizgen", "MERFISH", "Merscope" |
| AtlasXOmics | technology_docs/atlasxomics.md | "AtlasXOmics", "AtlasX", "spatial ATAC" |
| Takara | technology_docs/takara_workflow.md | "Takara", "Seeker", "Trekker" |
| **Latch APIs** | | |
| Remote files | latch_api_docs/lpath.md | "latch://", LPath, remote files |
| Widgets | latch_api_docs/plots_docs/widget-types.mdx | w_select, w_ldata_picker, user input |
| Plotting | latch_api_docs/plots_docs/custom-plots.mdx | matplotlib, plotly, w_plot |
| Signals | latch_api_docs/plots_docs/reactivity.mdx | Signal, reactivity, cross-cell dependencies |
