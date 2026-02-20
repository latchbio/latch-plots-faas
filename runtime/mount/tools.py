from __future__ import annotations

import json
from typing import Any

from claude_agent_sdk import create_sdk_mcp_server, tool
from lplots import _inject


VALID_NEXT_STATUSES = {
    "executing",
    "fixing",
    "thinking",
    "awaiting_user_response",
    "awaiting_cell_execution",
    "awaiting_user_widget_input",
    "done",
}

MCP_SERVER_NAME = "plots-agent-tools"

MCP_TOOL_NAMES = [
    "create_cell",
    "create_markdown_cell",
    "edit_cell",
    "run_cell",
    "update_plan",
    "submit_response",
    "delete_cell",
    "stop_cell",
    "delete_all_cells",
    "rename_notebook",
    "create_tab",
    "rename_tab",
    "restore_checkpoint",
    "execute_code",
    "get_global_info",
    "capture_widget_image",
    "set_widget",
    "h5_filter_by",
    "h5_color_by",
    "h5_refresh",
    "h5_set_selected_obsm_key",
    "h5_set_background_image",
    "h5_open_image_aligner",
    "h5_autoscale",
    "h5_zoom",
    "h5_set_background_image_visibility",
    "h5_add_selected_cells_to_categorical_obs",
    "h5_set_marker_opacity",
    "h5_manage_obs",
    "redeem_package",
    "smart_ui_spotlight",
]
MCP_ALLOWED_TOOL_NAMES = [f"mcp__{MCP_SERVER_NAME}__{name}" for name in MCP_TOOL_NAMES]


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


@tool(
    "delete_cell",
    "Remove a cell from the notebook.",
    {
        "type": "object",
        "properties": {
            "cell_id": {"type": "string", "description": "ID of the cell to delete"},
            "title": {"type": "string", "description": "Name of the cell to delete"},
            "action_summary": {
                "type": "string",
                "description": "Summary of the purpose of the delete.",
            },
        },
        "required": ["cell_id", "title", "action_summary"],
    },
)
async def delete_cell(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    cell_id = args["cell_id"]
    title = args["title"]
    action_summary = args["action_summary"]

    print(f"[tool] delete_cell id={cell_id}")
    result = await h.atomic_operation("delete_cell", {"cell_id": cell_id})
    if result.get("status") == "success":
        remaining = result.get("remaining_cells", [])
        cell_count = result.get("cell_count", 0)
        if remaining:
            cell_list = ", ".join([f"{c['index']}: {c['cell_type']}" for c in remaining[:5]])
            if len(remaining) > 5:
                cell_list += f", ... ({len(remaining) - 5} more)"
            msg = f"Cell {cell_id} deleted. {cell_count} cells remain: [{cell_list}]"
        else:
            msg = f"Cell {cell_id} deleted. No cells remain in notebook."
        print(f"[tool] delete_cell -> {msg}")
        return ok({
            "tool_name": "delete_cell",
            "summary": msg,
            "cell_id": cell_id,
            "cell_name": title,
            "message": action_summary,
            "success": True,
        })

    return ok({
        "tool_name": "delete_cell",
        "summary": f"Failed to delete cell: {result.get('error', 'Unknown error')}",
        "success": False,
    })


@tool(
    "stop_cell",
    "Stop execution of a specific cell.",
    {
        "type": "object",
        "properties": {
            "cell_id": {"type": "string", "description": "ID of the cell to stop"},
            "title": {"type": "string", "description": "Name of the cell to stop"},
            "action_summary": {
                "type": "string",
                "description": "Summary of the purpose of the stop.",
            },
        },
        "required": ["cell_id", "title", "action_summary"],
    },
)
async def stop_cell(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    cell_id = args["cell_id"]
    title = args["title"]
    action_summary = args["action_summary"]

    result = await h.atomic_operation("stop_cell", {"cell_id": cell_id})
    if result.get("status") == "success":
        h.executing_cells.discard(str(cell_id))
        return ok({
            "tool_name": "stop_cell",
            "summary": f"Stopped cell {cell_id}",
            "cell_id": cell_id,
            "cell_name": title,
            "message": action_summary,
            "success": True,
        })

    return ok({
        "tool_name": "stop_cell",
        "summary": f"Failed to stop cell {cell_id}: {result.get('error', 'Unknown error')}",
        "success": False,
    })


@tool(
    "delete_all_cells",
    "Delete all cells in the notebook efficiently.",
    {"type": "object", "properties": {}},
)
async def delete_all_cells(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    _ = args

    context_result = await h.atomic_operation("get_context", {})
    if context_result.get("status") != "success":
        return ok({
            "tool_name": "delete_all_cells",
            "summary": f"Failed to delete cells: {context_result.get('error', 'Unknown error')}",
            "success": False,
        })

    cells = context_result.get("context", {}).get("cells", [])
    deleted_count = 0
    for cell in reversed(cells):
        cell_id = cell.get("cell_id")
        if cell_id:
            result = await h.atomic_operation("delete_cell", {"cell_id": cell_id})
            if result.get("status") == "success":
                deleted_count += 1

    return ok({
        "tool_name": "delete_all_cells",
        "success": True,
        "summary": f"Deleted {deleted_count} cells from the notebook",
        "deleted_count": deleted_count,
    })


@tool(
    "rename_notebook",
    (
        "Rename the current plot notebook. Only call when the user explicitly asks "
        "to rename the notebook or the notebook is named Untitled Layout."
    ),
    {
        "type": "object",
        "properties": {
            "name": {
                "type": "string",
                "description": "New notebook name (<=5 words, Title Case)",
            }
        },
        "required": ["name"],
    },
)
async def rename_notebook(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    name = args["name"]
    print(f"[tool] rename_notebook name={name}")
    result = await h.atomic_operation("rename_notebook", {"name": name})
    if result.get("status") == "success":
        return ok({
            "tool_name": "rename_notebook",
            "success": True,
            "summary": f"Notebook renamed to '{name}'",
            "name": name,
        })

    return ok({
        "tool_name": "rename_notebook",
        "summary": f"Failed to rename notebook: {result.get('error', 'Unknown error')}",
        "success": False,
    })


@tool(
    "create_tab",
    (
        "Create a new tab marker cell at specified position to organize cells. "
        "This inserts a new cell and shifts all subsequent cell positions."
    ),
    {
        "type": "object",
        "properties": {
            "position": {
                "type": "integer",
                "description": "Position to insert the tab marker",
            },
            "display_name": {"type": "string", "description": "Name for the tab"},
        },
        "required": ["position", "display_name"],
    },
)
async def create_tab(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    position = args["position"]
    display_name = args["display_name"]
    if position < 0:
        return ok({
            "tool_name": "create_tab",
            "summary": "Error: Position must be non-negative",
            "success": False,
        })

    print(f'[tool] create_tab pos={position} name="{display_name}"')
    result = await h.atomic_operation(
        "create_tab", {"position": position, "display_name": display_name}
    )
    if result.get("status") == "success":
        tab_id = result.get("tab_id", "unknown")
        msg = f"Created tab at position {position} (ID: {tab_id}, Name: {display_name})"
        print(f"[tool] create_tab -> {msg}")
        return ok({
            "tool_name": "create_tab",
            "summary": msg,
            "tab_id": tab_id,
            "display_name": display_name,
            "position": position,
            "success": True,
        })

    return ok({
        "tool_name": "create_tab",
        "summary": f"Failed to create tab: {result.get('error', 'Unknown error')}",
        "success": False,
    })


@tool(
    "rename_tab",
    'Rename a tab. Use tab_id="DEFAULT" to rename the default tab.',
    {
        "type": "object",
        "properties": {
            "tab_id": {
                "type": "string",
                "description": (
                    'ID of the tab to rename. Use "DEFAULT" for the default tab, '
                    "or the TAB_ID from a Tab Marker."
                ),
            },
            "new_name": {"type": "string", "description": "New name for the tab"},
        },
        "required": ["tab_id", "new_name"],
    },
)
async def rename_tab(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    tab_id = args["tab_id"]
    new_name = args["new_name"]

    print(f'[tool] rename_tab tab_id={tab_id} new_name="{new_name}"')
    result = await h.atomic_operation("rename_tab", {"tab_id": tab_id, "new_name": new_name})
    if result.get("status") == "success":
        target = "Tab 1" if tab_id == "DEFAULT" else f"tab {tab_id}"
        msg = f"Renamed {target} to '{new_name}'"
        print(f"[tool] rename_tab -> {msg}")
        return ok({
            "tool_name": "rename_tab",
            "summary": msg,
            "tab_id": tab_id,
            "new_name": new_name,
            "success": True,
        })

    return ok({
        "tool_name": "rename_tab",
        "summary": f"Failed to rename tab: {result.get('error', 'Unknown error')}",
        "success": False,
    })


@tool(
    "restore_checkpoint",
    (
        "Restore this notebook to one of its own previously-created checkpoints "
        "(template_version_id)."
    ),
    {
        "type": "object",
        "properties": {
            "template_version_id": {
                "type": "string",
                "description": "Template version ID to restore",
            }
        },
        "required": ["template_version_id"],
    },
)
async def restore_checkpoint(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    template_version_id = args.get("template_version_id")
    print(f"[tool] restore_checkpoint template_version_id={template_version_id}")
    result = await h.atomic_operation(
        "restore_checkpoint", {"template_version_id": template_version_id}
    )
    if result.get("status") == "success":
        return ok({
            "tool_name": "restore_checkpoint",
            "success": True,
            "summary": f"Restored checkpoint to template version {template_version_id}",
            "template_version_id": template_version_id,
        })

    return ok({
        "tool_name": "restore_checkpoint",
        "summary": f"Failed to restore checkpoint: {result.get('error', 'Unknown error')}",
        "success": False,
    })


@tool(
    "execute_code",
    (
        "Execute arbitrary Python code in the notebook kernel and return the result, "
        "stdout, stderr, and any exceptions."
    ),
    {
        "type": "object",
        "properties": {"code": {"type": "string", "description": "Python code to execute"}},
        "required": ["code"],
    },
)
async def execute_code(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    code = args.get("code")
    if code is None:
        return ok({
            "tool_name": "execute_code",
            "success": False,
            "summary": "No code provided",
        })

    print(f"[tool] execute_code: {str(code)[:50]}...")
    result = await h.atomic_operation("execute_code", {"code": code})
    return ok({
        "tool_name": "execute_code",
        "success": result.get("status") == "success",
        "summary": (
            "Code executed"
            if result.get("status") == "success"
            else f"Failed to execute code: {result.get('error', 'Unknown error')}"
        ),
        "code": code,
        "stdout": result.get("stdout"),
        "stderr": result.get("stderr"),
        "exception": result.get("exception"),
    })


@tool(
    "get_global_info",
    (
        "Get rich information about a specific global variable including type, shape, "
        "columns, dtypes, and related metadata."
    ),
    {
        "type": "object",
        "properties": {
            "key": {
                "type": "string",
                "description": "Name of the global variable to inspect",
            }
        },
        "required": ["key"],
    },
)
async def get_global_info(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    key = args.get("key")
    if key is None:
        return ok({
            "tool_name": "get_global_info",
            "success": False,
            "summary": "No key provided",
        })

    print(f"[tool] get_global_info: {key}")
    result = await h.atomic_operation("get_global_info", {"key": key})
    if result.get("status") == "success":
        return ok({
            "tool_name": "get_global_info",
            "success": True,
            "summary": f"Retrieved info for global '{key}'",
            "key": key,
            "info": result.get("info", {}),
        })

    return ok({
        "tool_name": "get_global_info",
        "success": False,
        "summary": f"Failed to get global info: {result.get('error', 'Unknown error')}",
    })


@tool(
    "capture_widget_image",
    (
        "Capture a visual screenshot of an h5/AnnData or plot widget and return "
        "a base64-encoded PNG plus metadata."
    ),
    {
        "type": "object",
        "properties": {
            "widget_key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            }
        },
        "required": ["widget_key"],
    },
)
async def capture_widget_image(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    if widget_key is None:
        return ok({
            "tool_name": "capture_widget_image",
            "success": False,
            "summary": "No widget_key provided",
        })

    print(f"[tool] capture_widget_image: {widget_key}")
    result = await h.atomic_operation("capture_widget_image", {"widget_key": widget_key})
    if result.get("status") == "success":
        return ok({
            "tool_name": "capture_widget_image",
            "success": True,
            "summary": f"Captured image from {result.get('widget_type', 'unknown')} widget '{widget_key}'",
            "widget_key": widget_key,
            "widget_type": result.get("widget_type", "unknown"),
            "image": result.get("image"),
            "metadata": result.get("metadata", {}),
        })

    return ok({
        "tool_name": "capture_widget_image",
        "success": False,
        "summary": f"Failed to capture widget image: {result.get('error', 'Unknown error')}",
    })


@tool(
    "set_widget",
    "Set a single widget value by widget key.",
    {
        "type": "object",
        "properties": {
            "key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            },
            "value": {"description": "JSON-serializable value"},
            "action_summary": {
                "type": "string",
                "description": "Summary of the purpose of the set_widget.",
            },
            "label": {"type": "string", "description": "Label of the widget to set"},
        },
        "required": ["key", "value", "action_summary", "label"],
    },
)
async def set_widget(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    key = args.get("key")
    action_summary = args.get("action_summary")
    label = args.get("label")
    if not key:
        return ok({
            "tool_name": "set_widget",
            "summary": "Failed to set widget: Widget key is required",
            "success": False,
        })
    value = args.get("value")
    if value is None:
        return ok({
            "tool_name": "set_widget",
            "summary": "Failed to set widget: Widget value is required",
            "success": False,
        })

    print(f"[tool] set_widget key={key} value={value!r}")
    result = await h.atomic_operation("set_widget", {"key": key, "value": json.dumps(value)})
    if result.get("status") == "success":
        return ok({
            "tool_name": "set_widget",
            "success": True,
            "summary": f"Updated widget value for: {key}",
            "key": key,
            "label": label,
            "value": value,
            "message": action_summary,
        })

    return ok({
        "tool_name": "set_widget",
        "summary": f"Failed to update widget value: {result.get('error', 'Unknown error')}",
        "success": False,
    })


@tool(
    "h5_filter_by",
    (
        "Set filters for an h5/AnnData widget. Pass the complete array of filters "
        "to apply."
    ),
    {
        "type": "object",
        "properties": {
            "widget_key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            },
            "label": {"type": "string", "description": "Label of the widget"},
            "filters": {"type": "array", "description": "Complete filter array"},
        },
        "required": ["widget_key", "filters"],
    },
)
async def h5_filter_by(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    filters = args.get("filters")
    if isinstance(filters, str):
        try:
            filters = json.loads(filters)
        except json.JSONDecodeError:
            return ok({
                "tool_name": "h5_filter_by",
                "success": False,
                "summary": f"filters is invalid JSON: {filters!r}",
            })

    print(f"[tool] h5_filter_by widget_key={widget_key} filters={filters}")
    result = await h.atomic_operation("h5_filter_by", {"widget_key": widget_key, "filters": filters})
    if result.get("status") == "success":
        return ok({
            "tool_name": "h5_filter_by",
            "success": True,
            "label": args.get("label"),
            "summary": f"Applied filters to h5 widget: {filters}",
            "widget_key": widget_key,
            "filters": filters,
        })

    return ok({
        "tool_name": "h5_filter_by",
        "success": False,
        "summary": f"Failed to apply filters to h5 widget: {result.get('error', 'Unknown error')}",
    })


@tool(
    "h5_color_by",
    "Set an h5/AnnData widget to color by an observation or variable.",
    {
        "type": "object",
        "properties": {
            "widget_key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            },
            "label": {"type": "string", "description": "Label of the widget"},
            "color_by": {
                "description": "Color configuration (null, obs object, or var object)"
            },
        },
        "required": ["widget_key", "color_by"],
    },
)
async def h5_color_by(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    color_by = args.get("color_by")
    if isinstance(color_by, str):
        try:
            color_by = json.loads(color_by)
        except json.JSONDecodeError:
            return ok({
                "tool_name": "h5_color_by",
                "success": False,
                "summary": f"color_by is invalid JSON: {color_by!r}",
            })

    print(f"[tool] h5_color_by widget_key={widget_key} color_by={color_by}")
    result = await h.atomic_operation("h5_color_by", {"widget_key": widget_key, "color_by": color_by})
    if result.get("status") == "success":
        return ok({
            "tool_name": "h5_color_by",
            "success": True,
            "label": args.get("label"),
            "summary": f"Set h5 widget coloring: {color_by}",
            "widget_key": widget_key,
            "color_by": color_by,
        })

    return ok({
        "tool_name": "h5_color_by",
        "success": False,
        "summary": f"Failed to set h5 widget coloring: {result.get('error', 'Unknown error')}",
    })


@tool(
    "h5_refresh",
    "Refresh an h5/AnnData widget after state changes.",
    {
        "type": "object",
        "properties": {
            "widget_key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            },
            "label": {"type": "string", "description": "Label of the widget"},
        },
        "required": ["widget_key"],
    },
)
async def h5_refresh(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    print(f"[tool] h5_refresh widget_key={widget_key}")
    result = await h.atomic_operation("h5_refresh", {"widget_key": widget_key})
    if result.get("status") == "success":
        return ok({
            "tool_name": "h5_refresh",
            "success": True,
            "label": args.get("label"),
            "summary": f"Refreshed h5 widget {widget_key}",
            "widget_key": widget_key,
        })
    return ok({
        "tool_name": "h5_refresh",
        "success": False,
        "summary": f"Failed to refresh h5 widget: {result.get('error', 'Unknown error')}",
    })


@tool(
    "h5_set_selected_obsm_key",
    "Set the selected obsm key for an h5/AnnData widget.",
    {
        "type": "object",
        "properties": {
            "widget_key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            },
            "label": {"type": "string", "description": "Label of the widget"},
            "obsm_key": {"type": "string", "description": "obsm key to use"},
        },
        "required": ["widget_key", "obsm_key"],
    },
)
async def h5_set_selected_obsm_key(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    obsm_key = args.get("obsm_key")
    print(f"[tool] h5_set_selected_obsm_key widget_key={widget_key} obsm_key={obsm_key}")
    result = await h.atomic_operation(
        "h5_set_selected_obsm_key", {"widget_key": widget_key, "obsm_key": obsm_key}
    )
    if result.get("status") == "success":
        return ok({
            "tool_name": "h5_set_selected_obsm_key",
            "success": True,
            "label": args.get("label"),
            "summary": f"Set h5 widget to use obsm key {obsm_key}",
            "widget_key": widget_key,
            "obsm_key": obsm_key,
        })
    return ok({
        "tool_name": "h5_set_selected_obsm_key",
        "success": False,
        "summary": f"Failed to set h5 widget obsm key: {result.get('error', 'Unknown error')}",
    })


@tool(
    "h5_set_background_image",
    "Set a background image for an h5/AnnData widget from a user-attached file.",
    {
        "type": "object",
        "properties": {
            "label": {"type": "string", "description": "Label of the widget"},
            "widget_key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            },
            "node_id": {
                "type": "string",
                "description": "LData node ID of the image file to use as background",
            },
        },
        "required": ["widget_key", "node_id"],
    },
)
async def h5_set_background_image(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    node_id = args.get("node_id")
    print(f"[tool] h5_set_background_image widget_key={widget_key} node_id={node_id}")
    result = await h.atomic_operation(
        "h5_set_background_image", {"widget_key": widget_key, "node_id": node_id}
    )
    if result.get("status") == "success":
        return ok({
            "tool_name": "h5_set_background_image",
            "success": True,
            "label": args.get("label"),
            "summary": f"Set background image for h5 widget using {node_id}",
            "widget_key": widget_key,
            "node_id": node_id,
        })
    return ok({
        "tool_name": "h5_set_background_image",
        "success": False,
        "summary": f"Failed to set background image: {result.get('error', 'Unknown error')}",
    })


@tool(
    "h5_open_image_aligner",
    "Open the image alignment modal for a background image in an h5/AnnData widget.",
    {
        "type": "object",
        "properties": {
            "label": {"type": "string", "description": "Label of the widget"},
            "widget_key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            },
            "background_image_id": {
                "type": "string",
                "description": "ID of the background image to align",
            },
        },
        "required": ["widget_key", "background_image_id"],
    },
)
async def h5_open_image_aligner(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    background_image_id = args.get("background_image_id")
    print(
        "[tool] h5_open_image_aligner "
        f"widget_key={widget_key} background_image_id={background_image_id}"
    )
    result = await h.atomic_operation(
        "h5_open_image_aligner",
        {"widget_key": widget_key, "background_image_id": background_image_id},
    )
    if result.get("status") == "success":
        return ok({
            "tool_name": "h5_open_image_aligner",
            "success": True,
            "label": args.get("label"),
            "summary": f"Opened image aligner for background image {background_image_id}",
            "widget_key": widget_key,
            "background_image_id": background_image_id,
        })
    return ok({
        "tool_name": "h5_open_image_aligner",
        "success": False,
        "summary": f"Failed to open image aligner: {result.get('error', 'Unknown error')}",
    })


@tool(
    "h5_autoscale",
    "Reset the plotted view in an h5/AnnData widget to autoscaled data bounds.",
    {
        "type": "object",
        "properties": {
            "label": {"type": "string", "description": "Label of the widget"},
            "widget_key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            },
        },
        "required": ["widget_key"],
    },
)
async def h5_autoscale(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    print(f"[tool] h5_autoscale widget_key={widget_key}")
    result = await h.atomic_operation("h5_autoscale", {"widget_key": widget_key})
    if result.get("status") == "success":
        return ok({
            "tool_name": "h5_autoscale",
            "success": True,
            "label": args.get("label"),
            "summary": f"Autoscaled h5 widget {widget_key} to data bounds",
            "widget_key": widget_key,
        })
    return ok({
        "tool_name": "h5_autoscale",
        "success": False,
        "summary": f"Failed to autoscale h5 widget: {result.get('error', 'Unknown error')}",
    })


@tool(
    "h5_zoom",
    "Zoom the Plotly view in an h5/AnnData widget in or out.",
    {
        "type": "object",
        "properties": {
            "widget_key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            },
            "label": {"type": "string", "description": "Label of the widget"},
            "direction": {
                "type": "string",
                "enum": ["in", "out"],
                "description": "Zoom direction",
            },
            "percentage": {
                "type": "number",
                "minimum": 0,
                "description": "Optional percentage change",
            },
        },
        "required": ["widget_key", "direction"],
    },
)
async def h5_zoom(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    direction = args.get("direction")
    percentage = args.get("percentage")

    print(f"[tool] h5_zoom widget_key={widget_key} direction={direction} percentage={percentage}")
    params: dict[str, object] = {"widget_key": widget_key, "direction": direction}
    if percentage is not None:
        params["percentage"] = percentage

    result = await h.atomic_operation("h5_zoom", params)
    if result.get("status") == "success":
        zoom_desc = f"zoom {direction}"
        if percentage is not None:
            zoom_desc += f" by {percentage}%"
        return ok({
            "tool_name": "h5_zoom",
            "success": True,
            "label": args.get("label"),
            "summary": f"Applied {zoom_desc} to h5 widget {widget_key}",
            "widget_key": widget_key,
            "direction": direction,
            "percentage": percentage,
        })
    return ok({
        "tool_name": "h5_zoom",
        "success": False,
        "summary": f"Failed to zoom h5 widget: {result.get('error', 'Unknown error')}",
    })


@tool(
    "h5_set_background_image_visibility",
    "Show or hide a specific background image in an h5/AnnData widget.",
    {
        "type": "object",
        "properties": {
            "label": {"type": "string", "description": "Label of the widget"},
            "widget_key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            },
            "background_image_id": {
                "type": "string",
                "description": "ID of the background image to show or hide",
            },
            "hidden": {
                "type": "boolean",
                "description": "Whether to hide (true) or show (false)",
            },
        },
        "required": ["widget_key", "background_image_id", "hidden"],
    },
)
async def h5_set_background_image_visibility(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    background_image_id = args.get("background_image_id")
    hidden = args.get("hidden")

    print(
        "[tool] h5_set_background_image_visibility "
        f"widget_key={widget_key} background_image_id={background_image_id} hidden={hidden}"
    )
    result = await h.atomic_operation(
        "h5_set_background_image_visibility",
        {
            "widget_key": widget_key,
            "background_image_id": background_image_id,
            "hidden": hidden,
        },
    )
    if result.get("status") == "success":
        visibility_action = "hidden" if hidden else "shown"
        return ok({
            "tool_name": "h5_set_background_image_visibility",
            "success": True,
            "label": args.get("label"),
            "summary": f"Background image {background_image_id} {visibility_action}",
            "widget_key": widget_key,
            "background_image_id": background_image_id,
            "hidden": hidden,
        })
    return ok({
        "tool_name": "h5_set_background_image_visibility",
        "success": False,
        "summary": f"Failed to set background image visibility: {result.get('error', 'Unknown error')}",
    })


@tool(
    "h5_add_selected_cells_to_categorical_obs",
    "Assign selected cells to a category in a categorical observation key.",
    {
        "type": "object",
        "properties": {
            "label": {"type": "string", "description": "Label of the widget"},
            "widget_key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            },
            "obs_key": {"type": "string", "description": "Categorical observation key"},
            "category": {
                "type": "string",
                "description": "Category name to assign to selected cells",
            },
        },
        "required": ["widget_key", "obs_key", "category"],
    },
)
async def h5_add_selected_cells_to_categorical_obs(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    obs_key = args.get("obs_key")
    category = args.get("category")

    print(
        "[tool] h5_add_selected_cells_to_categorical_obs "
        f"widget_key={widget_key} obs_key={obs_key} category={category}"
    )
    result = await h.atomic_operation(
        "h5_add_selected_cells_to_categorical_obs",
        {"widget_key": widget_key, "obs_key": obs_key, "category": category},
    )
    if result.get("status") == "success":
        return ok({
            "tool_name": "h5_add_selected_cells_to_categorical_obs",
            "success": True,
            "label": args.get("label"),
            "summary": f"Assigned selected cells to category '{category}' in observation key '{obs_key}'",
            "widget_key": widget_key,
            "obs_key": obs_key,
            "category": category,
        })
    return ok({
        "tool_name": "h5_add_selected_cells_to_categorical_obs",
        "success": False,
        "summary": f"Failed to assign selected cells to category: {result.get('error', 'Unknown error')}",
    })


@tool(
    "h5_set_marker_opacity",
    "Set marker opacity for all cell markers in an h5/AnnData widget.",
    {
        "type": "object",
        "properties": {
            "label": {"type": "string", "description": "Label of the widget"},
            "widget_key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            },
            "opacity": {
                "type": "number",
                "description": "Opacity value between 0.1 and 0.9",
            },
        },
        "required": ["widget_key", "opacity"],
    },
)
async def h5_set_marker_opacity(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    opacity = args.get("opacity")

    print(f"[tool] h5_set_marker_opacity widget_key={widget_key} opacity={opacity}")
    result = await h.atomic_operation(
        "h5_set_marker_opacity", {"widget_key": widget_key, "opacity": opacity}
    )
    if result.get("status") == "success":
        return ok({
            "tool_name": "h5_set_marker_opacity",
            "success": True,
            "label": args.get("label"),
            "summary": f"Set marker opacity to {opacity}",
            "widget_key": widget_key,
            "opacity": opacity,
        })
    return ok({
        "tool_name": "h5_set_marker_opacity",
        "success": False,
        "summary": f"Failed to set marker opacity: {result.get('error', 'Unknown error')}",
    })


@tool(
    "h5_manage_obs",
    "Create or delete an observation column in an h5/AnnData widget.",
    {
        "type": "object",
        "properties": {
            "widget_key": {
                "type": "string",
                "description": "Full widget key in format <tf_id>/<widget_id>",
            },
            "label": {"type": "string", "description": "Label of the widget"},
            "obs_key": {
                "type": "string",
                "description": "Observation column name to create or delete",
            },
            "operation": {
                "type": "string",
                "enum": ["add", "remove"],
                "description": "Whether to add or remove the observation column",
            },
            "obs_type": {
                "type": "string",
                "enum": ["category", "bool", "int64", "float64"],
                "description": "Observation type when operation is add",
            },
        },
        "required": ["widget_key", "obs_key", "operation"],
    },
)
async def h5_manage_obs(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    obs_key = args.get("obs_key")
    operation = args.get("operation")
    obs_type = args.get("obs_type", "category")

    print(
        "[tool] h5_manage_obs "
        f"widget_key={widget_key} obs_key={obs_key} operation={operation} obs_type={obs_type}"
    )
    params: dict[str, object] = {
        "widget_key": widget_key,
        "obs_key": obs_key,
        "operation": operation,
    }
    if operation == "add":
        params["obs_type"] = obs_type

    result = await h.atomic_operation("h5_manage_obs", params)
    if result.get("status") == "success":
        if operation == "add":
            return ok({
                "tool_name": "h5_manage_obs",
                "success": True,
                "label": args.get("label"),
                "summary": f"Created observation column '{obs_key}' with type '{obs_type}'",
                "widget_key": widget_key,
                "obs_key": obs_key,
                "operation": operation,
                "obs_type": obs_type,
            })
        return ok({
            "tool_name": "h5_manage_obs",
            "success": True,
            "label": args.get("label"),
            "summary": f"Deleted observation column '{obs_key}'",
            "widget_key": widget_key,
            "obs_key": obs_key,
            "operation": operation,
        })

    return ok({
        "tool_name": "h5_manage_obs",
        "success": False,
        "summary": f"Failed to {operation} observation column: {result.get('error', 'Unknown error')}",
    })


@tool(
    "redeem_package",
    "Redeem a package so the workspace gains access to technology-specific assets.",
    {
        "type": "object",
        "properties": {
            "package_code": {"type": "string", "description": "Multi-use package invite code"},
            "package_version_id": {"type": "string", "description": "Package version ID"},
            "redemption_reason": {
                "type": "string",
                "description": "Reason for redeeming the package",
            },
        },
        "required": ["package_code", "package_version_id", "redemption_reason"],
    },
)
async def redeem_package(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    package_code = args.get("package_code")
    package_version_id = args.get("package_version_id")
    redemption_reason = args.get("redemption_reason")
    if package_code is None or package_version_id is None or redemption_reason is None:
        return ok({
            "tool_name": "redeem_package",
            "success": False,
            "summary": (
                "Package code, package version ID and redemption reason are required "
                "to redeem a package"
            ),
        })

    package_code_str = str(package_code)
    package_version_id_str = str(package_version_id)
    print("[tool] redeem_package")
    result = await h.atomic_operation(
        "redeem_package",
        {"package_code": package_code_str, "package_version_id": package_version_id_str},
    )
    if result.get("status") == "success":
        return ok({
            "tool_name": "redeem_package",
            "success": True,
            "summary": f"Redeemed package {package_code_str} (version {package_version_id_str})",
            "package_code": package_code_str,
            "package_version_id": package_version_id_str,
            "redemption_reason": redemption_reason,
        })

    return ok({
        "tool_name": "redeem_package",
        "success": False,
        "summary": f"Failed to redeem package: {result.get('error', 'Unknown error')}",
    })


@tool(
    "smart_ui_spotlight",
    "Highlight a UI element to guide the user's attention.",
    {
        "type": "object",
        "properties": {
            "keyword": {
                "type": "string",
                "enum": ["lasso_select", "file_upload", "widget_input"],
                "description": "The UI element to highlight",
            },
            "widget_key": {
                "type": "string",
                "description": "Optional full widget key for widget-related highlights",
            },
            "widget_label": {
                "type": "string",
                "description": "Optional label of the widget to highlight",
            },
        },
        "required": ["keyword"],
    },
)
async def smart_ui_spotlight(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    keyword = args.get("keyword")
    widget_key = args.get("widget_key")
    widget_label = args.get("widget_label")

    print(f"[tool] smart_ui_spotlight keyword={keyword}, widget_key={widget_key}")
    params: dict[str, object] = {"keyword": keyword}
    if widget_key is not None:
        params["widget_key"] = widget_key

    result = await h.atomic_operation("smart_ui_spotlight", params)
    if result.get("status") == "success":
        payload: dict[str, object] = {
            "tool_name": "smart_ui_spotlight",
            "success": True,
            "summary": f"Highlighted UI element: {keyword}",
            "keyword": keyword,
        }
        if widget_key is not None:
            payload["widget_key"] = widget_key
        if widget_label is not None:
            payload["widget_label"] = widget_label
        return ok(payload)

    return ok({
        "tool_name": "smart_ui_spotlight",
        "success": False,
        "summary": f"Failed to highlight UI element: {result.get('error', 'Unknown error')}",
    })


@tool(
    "update_plan",
    "Update the in-memory execution plan shown in the agent UI.",
    {
        "type": "object",
        "properties": {
            "plan": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "id": {"type": "string"},
                        "description": {"type": "string"},
                        "status": {
                            "type": "string",
                            "enum": ["todo", "in_progress", "done", "cancelled"],
                        },
                    },
                    "required": ["id", "description", "status"],
                },
            },
            "plan_diff": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "id": {"type": "string"},
                        "action": {
                            "type": "string",
                            "enum": ["add", "update", "complete", "remove"],
                        },
                        "description": {"type": "string"},
                    },
                    "required": ["id", "action"],
                },
            },
            "plan_update_overview": {"type": "string"},
        },
        "required": ["plan", "plan_diff", "plan_update_overview"],
    },
)
async def update_plan(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    try:
        plan_items = args.get("plan", [])
        plan_diff = args.get("plan_diff", [])
        plan_update_overview = args.get("plan_update_overview", "")

        if not isinstance(plan_items, list):
            return ok({
                "tool_name": "update_plan",
                "success": False,
                "summary": "Invalid plan: expected list",
            })
        if not isinstance(plan_diff, list):
            return ok({
                "tool_name": "update_plan",
                "success": False,
                "summary": "Invalid plan_diff: expected list",
            })

        h.current_plan = {"steps": plan_items}

        print(f"[tool] update_plan: {plan_update_overview}")
        for item in plan_items:
            if isinstance(item, dict):
                print(
                    f"  [{item.get('status')}] {item.get('id')}: {item.get('description')}"
                )
        for diff in plan_diff:
            if isinstance(diff, dict):
                print(f"  diff: [{diff.get('action')}] {diff.get('id')}")

        summary = (
            plan_update_overview if isinstance(plan_update_overview, str) and plan_update_overview.strip() != "" else "Plan updated"
        )
        return ok({
            "tool_name": "update_plan",
            "success": True,
            "summary": summary,
        })
    except Exception as e:
        print(f"[tool] update_plan error: {e!s}")
        return ok({
            "tool_name": "update_plan",
            "success": False,
            "summary": f"Error updating plan: {e!s}",
        })


@tool(
    "submit_response",
    (
        "Submit user-facing structured response and set next_status. "
        "At least one of summary or questions is required."
    ),
    {
        "type": "object",
        "properties": {
            "summary": {"type": "string"},
            "questions": {"type": "string"},
            "next_status": {"type": "string", "enum": sorted(VALID_NEXT_STATUSES)},
            "expected_widgets": {
                "type": "array",
                "items": {"type": "string"},
            },
            "continue": {"type": "boolean", "default": False},
        },
        "required": ["next_status"],
    },
)
async def submit_response(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    try:
        summary = args.get("summary")
        if summary is not None and not isinstance(summary, str):
            summary = None

        questions = args.get("questions")
        if questions is not None and not isinstance(questions, str):
            questions = None

        if (summary is None or summary.strip() == "") and (
            questions is None or questions.strip() == ""
        ):
            return ok({
                "tool_name": "submit_response",
                "success": False,
                "summary": "Please provide at least one of summary or questions.",
            })

        next_status = args.get("next_status")
        if not isinstance(next_status, str) or next_status not in VALID_NEXT_STATUSES:
            return ok({
                "tool_name": "submit_response",
                "success": False,
                "summary": "Please provide a valid next_status.",
            })

        should_continue = bool(args.get("continue", False))

        expected_widgets_arg = args.get("expected_widgets", [])
        expected_widgets: list[str] = []
        if isinstance(expected_widgets_arg, list):
            expected_widgets = [
                str(widget_key)
                for widget_key in expected_widgets_arg
                if isinstance(widget_key, str) and widget_key.strip() != ""
            ]

        print("[tool] submit_response called with:")
        print(f"  - next_status: {next_status}")
        print(f"  - summary: {summary}")
        print(f"  - questions: {questions}")
        print(f"  - continue: {should_continue}")
        print(f"  - expected_widgets: {expected_widgets}")
        if should_continue:
            print(
                "[tool] submit_response: continue=true retained for compatibility; "
                "SDK runtime controls actual turn continuation."
            )

        if should_continue and h.executing_cells:
            h.should_auto_continue = False
            h.pending_auto_continue = True
        else:
            h.should_auto_continue = should_continue
            h.pending_auto_continue = False

        terminal_statuses = {"done", "awaiting_user_response"}
        h.pause_until_user_query = (
            not should_continue and next_status in terminal_statuses
        )

        h.current_status = next_status
        if next_status == "awaiting_user_widget_input":
            h.expected_widgets = {widget_key: None for widget_key in expected_widgets}
        else:
            h.expected_widgets.clear()

        return ok({
            "tool_name": "submit_response",
            "success": True,
            "summary": "Response submitted successfully",
            "next_status": next_status,
        })
    except Exception as e:
        print(f"[tool] submit_response error: {e!s}")
        return ok({
            "tool_name": "submit_response",
            "success": False,
            "summary": f"Error submitting response: {e!s}",
        })

# todo(rteqs): idk splitting them actually helps, i.e. one big mcp server or multiple small. maybe for subagents that we want to limit capability so context window stays small
core_tools = [
    create_cell,
    create_markdown_cell,
    edit_cell,
    delete_cell,
    run_cell,
    stop_cell,
    delete_all_cells,
    rename_notebook,
    create_tab,
    rename_tab,
    restore_checkpoint,
    execute_code,
    get_global_info,
    capture_widget_image,
    update_plan,
    submit_response,
]
widget_tools = [
    set_widget,
    h5_filter_by,
    h5_color_by,
    h5_refresh,
    h5_set_selected_obsm_key,
    h5_set_background_image,
    h5_open_image_aligner,
    h5_autoscale,
    h5_zoom,
    h5_set_background_image_visibility,
    h5_add_selected_cells_to_categorical_obs,
    h5_set_marker_opacity,
    h5_manage_obs,
]
misc_tools = [
    redeem_package,
    smart_ui_spotlight,
]

all_tools = core_tools + widget_tools + misc_tools


agent_tools_mcp = create_sdk_mcp_server(name=MCP_SERVER_NAME, tools=all_tools)
