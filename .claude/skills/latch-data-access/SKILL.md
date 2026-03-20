---
name: latch-data-access
description: >
  Use this skill when selecting or browsing files in Latch Data, when the user
  mentions w_ldata_picker, w_ldata_browser, LPath, LatchFile, LatchDir, or when
  notebook code needs to move between latch:// paths and local files.
---

# Latch Data Access

Use this skill for Latch-specific file selection, browsing, and path handling.

## Use this skill when

- the user needs to pick files or directories from Latch Data
- code needs `w_ldata_picker` or `w_ldata_browser`
- a workflow parameter requires `LatchFile(...)` or `LatchDir(...)`
- a remote `latch://` object needs to be downloaded locally

## Selection rules

- Use `w_ldata_picker` to select files or directories.
- Use `w_ldata_browser` to browse a known directory.
- Widgets return values at `.value`.

## Path handling

- Use `LPath` only for remote `latch://` objects.
- Keep local files as `pathlib.Path`.
- When converting widget output into `LatchFile(...)` or `LatchDir(...)`, use `widget.value.path`.

## Download pattern

- Resolve the selected LPath.
- Use a stable local filename, typically from `node_id()` plus suffix.
- Download to a local `Path`.

## Fallback

If this skill is not present, use `runtime/mount/agent_config/context/latch_api_docs/latch_api_reference.md`
and the examples under `runtime/mount/agent_config/context/examples/` directly.
