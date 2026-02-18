from __future__ import annotations

import json
from typing import Any

from claude_agent_sdk import create_sdk_mcp_server, tool
from lplots import _inject

def ok(d: dict) -> dict[str, Any]:
    return {"content": [{"type": "text", "text": json.dumps(d)}]}


def harness():
    h = _inject.agent
    if h is None:
        raise RuntimeError("Agent harness is not initialized")
    return h


@tool(
    "create_cell",
    "Create a new code cell at specified position. The cell will automatically run after creation.",
    {
        "type": "object",
        "properties": {
            "position": {
                "type": "integer",
                "description": "Position to insert the cell",
            },
            "code": {"type": "string", "description": "Python code for the cell"},
            "title": {"type": "string", "description": "Name for the cell"},
            "action_summary": {
                "type": "string",
                "description": "Summary of the purpose of the cell.",
            },
        },
        "required": ["position", "code", "title", "action_summary"],
    },
)
async def create_cell(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    position = args["position"]
    code = args["code"]
    title = args["title"]
    action_summary = args["action_summary"]

    if position < 0:
        return ok({
            "tool_name": "create_cell",
            "summary": "Error: Position must be non-negative",
            "success": False,
        })

    print(f'[tool] create_cell pos={position} title="{title}"')
    params = {
        "position": position,
        "cell_type": "code",
        "source": code,
        "title": title,
        "auto_run": True,
    }
    result = await h.atomic_operation("create_cell", params)

    if result.get("status") == "success":
        cell_id = result.get("cell_id", "unknown")
        tf_id = result.get("tf_id", "unknown")
        msg = f"Created cell at position {position} (cell_id: {cell_id}, tf_id: {tf_id}, title: {title})"
        print(f"[tool] create_cell -> {msg}")
        return ok({
            "tool_name": "create_cell",
            "summary": msg,
            "code": code,
            "cell_id": cell_id,
            "tf_id": tf_id,
            "cell_name": title,
            "position": position,
            "message": action_summary,
            "success": True,
        })
    return ok({
        "tool_name": "create_cell",
        "summary": f"Failed to create cell: {result.get('error', 'Unknown error')}",
        "success": False,
    })


@tool(
    "create_markdown_cell",
    "Create a new markdown cell at specified position.",
    {
        "type": "object",
        "properties": {
            "position": {
                "type": "integer",
                "description": "Position to insert the cell",
            },
            "code": {"type": "string", "description": "Markdown content"},
            "title": {
                "type": "string",
                "description": "Title of first header in the markdown cell",
            },
            "action_summary": {
                "type": "string",
                "description": "Summary of the purpose of the cell.",
            },
        },
        "required": ["position", "code", "title", "action_summary"],
    },
)
async def create_markdown_cell(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    position = args["position"]
    code = args["code"]
    title = args["title"]
    action_summary = args["action_summary"]

    if position < 0:
        return ok({"summary": "Error: Position must be non-negative", "success": False})

    print(f"[tool] create_markdown_cell pos={position}")
    params = {"position": position, "cell_type": "markdown", "source": code}
    result = await h.atomic_operation("create_markdown_cell", params)

    if result.get("status") == "success":
        cell_id = result.get("cell_id", "unknown")
        msg = f"Created markdown cell at position {position} (ID: {cell_id})"
        print(f"[tool] create_markdown_cell -> {msg}")
        return ok({
            "tool_name": "create_markdown_cell",
            "summary": msg,
            "code": code,
            "cell_id": cell_id,
            "cell_name": title,
            "position": position,
            "message": action_summary,
            "success": True,
        })
    return ok({
        "tool_name": "create_markdown_cell",
        "message": f"Failed to create cell: {result.get('error', 'Unknown error')}",
        "success": False,
    })


@tool(
    "edit_cell",
    "Replace the contents of an existing cell. The cell will automatically run after editing.",
    {
        "type": "object",
        "properties": {
            "cell_id": {"type": "string", "description": "ID of the cell to edit"},
            "new_code": {
                "type": "string",
                "description": "New code/content for the cell",
            },
            "title": {"type": "string", "description": "Name of the cell to edit"},
            "action_summary": {
                "type": "string",
                "description": "Summary of the purpose of the edit.",
            },
        },
        "required": ["cell_id", "new_code", "title"],
    },
)
async def edit_cell(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    cell_id = args["cell_id"]
    new_code = args["new_code"]
    title = args["title"]
    action_summary = args.get("action_summary", "")

    print(f"[tool] edit_cell id={cell_id}")
    original_code = ""
    for cell in h.latest_notebook_context.get("cells", []):
        if cell.get("cell_id") == cell_id:
            original_code = cell.get("source", "")
            break

    params = {"cell_id": cell_id, "source": new_code, "auto_run": True}
    result = await h.atomic_operation("edit_cell", params)

    if result.get("status") == "success":
        msg = f"Cell {cell_id} edited successfully"
        print(f"[tool] edit_cell -> {msg}")
        return ok({
            "tool_name": "edit_cell",
            "summary": msg,
            "code": new_code,
            "original_code": original_code,
            "cell_id": cell_id,
            "cell_name": title,
            "message": action_summary,
            "success": True,
        })
    return ok({
        "tool_name": "edit_cell",
        "summary": f"Failed to edit cell: {result.get('error', 'Unknown error')}",
        "success": False,
    })


@tool(
    "run_cell",
    "Run an existing code cell.",
    {
        "type": "object",
        "properties": {
            "cell_id": {"type": "string", "description": "ID of the cell to run"},
            "title": {"type": "string", "description": "Name of the cell to run"},
            "action_summary": {
                "type": "string",
                "description": "Summary of the purpose of the run.",
            },
        },
        "required": ["cell_id", "title", "action_summary"],
    },
)
async def run_cell(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    cell_id = args["cell_id"]
    title = args["title"]
    action_summary = args["action_summary"]

    print(f"[tool] run_cell id={cell_id} title={title!r}")
    await h.send({
        "type": "agent_action",
        "action": "run_cell",
        "params": {"cell_id": cell_id},
    })
    h.executing_cells.add(str(cell_id))

    return ok({
        "tool_name": "run_cell",
        "summary": f"Cell {cell_id} execution started",
        "cell_id": cell_id,
        "cell_name": title,
        "message": action_summary,
        "success": True,
    })

# todo(rteqs): every other tool in the current harness

# todo(rteqs): idk splitting them actually helps, i.e. one big mcp server or multiple small. maybe for subagents that we want to limit capability so context window stays small
notebook_tools = [
    create_cell,
    create_markdown_cell,
    edit_cell,
    run_cell,
    # delete_cell,
    # run_cell,
    # stop_cell,
    # delete_all_cells,
    # rename_notebook,
    # create_tab,
    # rename_tab,
    # restore_checkpoint,
    # get_global_info
    # execute_code
]
widget_tools = [
    # set_widget,
    # h5_filter_by,
    # h5_color_by,
    # h5_refresh,
    # h5_set_selected_obsm_key,
    # h5_set_background_image,
    # h5_open_image_aligner,
    # h5_autoscale,
    # h5_zoom,
    # h5_set_background_image_visibility,
    # h5_add_selected_cells_to_categorical_obs,
    # h5_set_marker_opacity,
    # h5_manage_obs,
    # capture_widget_image
]
misc_tools = [
    # update_plan,
    # submit_response,
    # redeem_package,
    # smart_ui_spotlight
]

all_tools = notebook_tools + widget_tools + misc_tools


notebook_mcp = create_sdk_mcp_server(name="notebook-tools", tools=notebook_tools)
# todo(rteqs): other mcp servers
