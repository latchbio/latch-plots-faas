---
name: latch-workflows
description: >
  Use this skill when launching or inspecting Latch workflows in Plots, when the
  user mentions w_workflow, workflow params, LatchFile, LatchDir, workflow outputs,
  or waiting for workflow completion.
---

# Latch Workflows

Use this skill for Latch-specific workflow execution mechanics.

## Use this skill when

- the plan requires `w_workflow`
- a technology doc references a workflow in `wf/`
- the user needs help constructing workflow params
- the user needs to validate `LatchFile(...)` or `LatchDir(...)` inputs
- the workflow has been launched and you need to wait for outputs

## Before launching

- Read the technology-specific workflow reference first if one is available.
- Build `params` exactly as the workflow reference specifies.
- When a widget returns an LPath-like value, use `.path` when constructing `LatchFile(...)` or `LatchDir(...)`.
- Show the complete `params` before launch using markdown or `w_text_output`.

## Required validation

Before calling `w_workflow`, verify:

- no `None` values
- no empty `LatchFile()` or `LatchDir()`
- all required parameters are present
- all Latch paths are valid strings

If any parameter is invalid, stop and fix the cell before proceeding.

## Launch pattern

- Always call `w_workflow(..., automatic=True)`.
- Always use a unique `key`.
- After launch, retrieve `execution = w.value`.
- If execution exists, `await execution.wait()`.
- Inspect outputs only after the workflow reaches a terminal state.

## Fallback

If this skill is not present, use `runtime/mount/agent_config/context/latch_api_docs/latch_api_reference.md`
and the technology repo's local `wf/` docs directly.
