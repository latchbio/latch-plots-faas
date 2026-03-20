---
name: latch-plots-ui
description: >
  Use this skill when working with Latch plot or output widgets in Plots, including
  w_h5, w_plot, w_table, w_text_output, viewer refresh behavior, or notebook UI
  rendering conventions.
---

# Latch Plots UI

Use this skill for Latch-specific rendering and notebook UI conventions.

## Use this skill when

- the plan needs `w_h5`, `w_plot`, `w_table`, or `w_text_output`
- the user asks how to render AnnData, tables, or plots in Plots
- a technology doc requires plot widgets or viewer output

## Output rules

- Use `w_text_output` for short status or validation messages.
- Use `w_table` for DataFrames.
- Use `w_plot` for Plotly or Matplotlib figures.
- Use `w_h5` for AnnData exploration and spatial / embedding inspection.

## `w_table`

- Pass a named DataFrame variable directly.
- Do not pass expressions like `df.head()` directly to `w_table`.

## `w_plot`

- Each plot must use its own unique named variable.
- Do not rely on `globals()` or dynamic variable names inside loops.
- For Scanpy dot or violin plots, convert the returned object into a figure before passing it to `w_plot`.

## `w_h5`

- Use `w_h5` for interactive AnnData exploration.
- Keep existing viewers in sync after updating `adata.obs` or `adata.obsm`.

## Fallback

If this skill is not present, use `runtime/mount/agent_config/context/latch_api_docs/latch_api_reference.md`
and the examples under `runtime/mount/agent_config/context/examples/` directly.
