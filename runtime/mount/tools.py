import json
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from claude_agent_sdk import (
    PermissionResultAllow,
    PermissionResultDeny,
    ToolPermissionContext,
    create_sdk_mcp_server,
    tool,
)
from latch_data_validation.data_validation import validate
from loro import LoroMap
from lplots import _inject
from notebook import get_notebook, parse_ts_container_id
from utils import auth_token_sdk, gql_query


@dataclass(frozen=True)
class RedeemPackageResult:
    resPackageId: str
    resPackageVersionId: str
    resRedeemedInAccountId: str


@dataclass(frozen=True)
class RedeemPackagePayload:
    clientMutationId: str | None
    result: RedeemPackageResult | None


@dataclass(frozen=True)
class RedeemPackageData:
    redeemPackage: RedeemPackagePayload | None


@dataclass(frozen=True)
class RedeemPackageRes:
    data: RedeemPackageData


@dataclass(frozen=True)
class PlotTransformInfoNode:
    id: str


@dataclass(frozen=True)
class PlotTransformInfoNodes:
    nodes: list[PlotTransformInfoNode]


@dataclass(frozen=True)
class ResolveTfIdData:
    plotTransformInfos: PlotTransformInfoNodes


@dataclass(frozen=True)
class ResolveTfIdRes:
    data: ResolveTfIdData


async def resolve_tf_id(notebook_id: str, loro_cell_id: str) -> str | None:
    gql_res = await gql_query(
        query="""
            query ResolveTfId($notebookId: BigInt!, $loroCellId: String!) {
                plotTransformInfos(
                    filter: {
                        parentNotebookId: { equalTo: $notebookId }
                        loroCellId: { equalTo: $loroCellId }
                        removalTime: { isNull: true }
                    }
                    first: 1
                ) {
                    nodes {
                        id
                    }
                }
            }
        """,
        variables={"notebookId": notebook_id, "loroCellId": loro_cell_id},
        auth=auth_token_sdk,
    )
    res = validate(gql_res, ResolveTfIdRes)
    nodes = res.data.plotTransformInfos.nodes
    if len(nodes) == 0:
        return None
    return nodes[0].id


if TYPE_CHECKING:
    from agent import AgentHarness

MCP_SERVER_NAME = "plots-agent-tools"

MCP_TOOL_NAMES = [
    "create_cell",
    "create_markdown_cell",
    "edit_cell",
    "run_cell",
    "update_plan",
    "delete_cell",
    "stop_cell",
    "delete_all_cells",
    "rename_notebook",
    "create_tab",
    "rename_tab",
    "restore_checkpoint",
    "execute_code",
    "get_global_info",
    "set_widget",
    "get_widget",
    "get_plot_image",
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
    "submit_response",
]
MCP_ALLOWED_TOOL_NAMES = [f"mcp__{MCP_SERVER_NAME}__{name}" for name in MCP_TOOL_NAMES]


def ok(d: dict) -> dict[str, Any]:
    return {"content": [{"type": "text", "text": json.dumps(d)}]}


def harness() -> "AgentHarness":
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
    try:
        notebook = await get_notebook()
        cell_id, tf_id = await notebook.create_code_cell(
            pos=position, code=code, display_name=title
        )
    except Exception as e:
        return ok({
            "tool_name": "create_cell",
            "summary": f"Failed to create cell: {e!s}",
            "success": False,
        })

    h.pending_cells.add(tf_id)
    await h.send({
        "type": "agent_action",
        "action": "run_cell",
        "params": {"cell_id": tf_id, "code": code, "parallel": False},
    })

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
    position = args["position"]
    code = args["code"]
    title = args["title"]
    action_summary = args["action_summary"]

    if position < 0:
        return ok({"summary": "Error: Position must be non-negative", "success": False})

    print(f"[tool] create_markdown_cell pos={position}")
    try:
        notebook = await get_notebook()
        cell_id = await notebook.create_markdown_cell(pos=position, content=code)
    except Exception as e:
        return ok({
            "tool_name": "create_markdown_cell",
            "message": f"Failed to create cell: {e!s}",
            "success": False,
        })

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
    action_summary = args["action_summary"]

    print(f"[tool] edit_cell id={cell_id}")
    try:
        notebook = await get_notebook()
        original_code = ""
        container = notebook.loro_doc.get_container(id=parse_ts_container_id(cell_id))
        if isinstance(container, LoroMap):
            deep = container.get_deep_value()
            if isinstance(deep, dict):
                original_code = str(deep.get("source", "") or "")

        tf_id = await resolve_tf_id(notebook.notebook_id, cell_id)
        if tf_id is None:
            return ok({
                "tool_name": "edit_cell",
                "summary": f"Failed to edit cell: no transform found for cell_id {cell_id}",
                "success": False,
            })

        await notebook.edit_cell(cell_id=cell_id, new_code=new_code)
    except Exception as e:
        return ok({
            "tool_name": "edit_cell",
            "summary": f"Failed to edit cell: {e!s}",
            "success": False,
        })

    h.pending_cells.add(tf_id)
    await h.send({
        "type": "agent_action",
        "action": "run_cell",
        "params": {"cell_id": tf_id, "code": new_code, "parallel": False},
    })

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


# todo(rteqs): when ClaudeCode support MCP tasks (https://modelcontextprotocol.io/seps/1686-tasks) we should make this a background task and handle tasks from SystemMessges
@tool(
    "run_cell",
    "Execute a specific cell.",
    {
        "type": "object",
        "properties": {
            "cell_id": {
                "type": "string",
                "description": "ID of the cell to run (note: loro_cell_id, not tf_id)",
            },
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

    print(f"[tool] run_cell id={cell_id}")
    try:
        notebook = await get_notebook()
        container = notebook.loro_doc.get_container(id=parse_ts_container_id(cell_id))
        if not isinstance(container, LoroMap):
            return ok({
                "tool_name": "run_cell",
                "summary": f"Failed to run cell: cell {cell_id} not found",
                "success": False,
            })
        deep = container.get_deep_value()
        source = ""
        if isinstance(deep, dict):
            source = str(deep.get("source", "") or "")

        tf_id = await resolve_tf_id(notebook.notebook_id, cell_id)
        if tf_id is None:
            return ok({
                "tool_name": "run_cell",
                "summary": f"Failed to run cell: no transform found for cell_id {cell_id}",
                "success": False,
            })
    except Exception as e:
        return ok({
            "tool_name": "run_cell",
            "summary": f"Failed to run cell: {e!s}",
            "success": False,
        })

    h.pending_cells.add(tf_id)
    await h.send({
        "type": "agent_action",
        "action": "run_cell",
        "params": {"cell_id": tf_id, "code": source, "parallel": False},
    })

    return ok({
        "tool_name": "run_cell",
        "summary": f"Cell {cell_id} execution started",
        "cell_id": cell_id,
        "tf_id": tf_id,
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
    cell_id = args["cell_id"]
    title = args["title"]
    action_summary = args["action_summary"]

    print(f"[tool] delete_cell id={cell_id}")
    try:
        notebook = await get_notebook()
        await notebook.delete_cell(cell_id=cell_id)
    except Exception as e:
        return ok({
            "tool_name": "delete_cell",
            "summary": f"Failed to delete cell: {e!s}",
            "success": False,
        })

    remaining: list[dict[str, Any]] = []
    for i, cid in enumerate(notebook.cells.to_vec()):
        container = notebook.loro_doc.get_container(id=cid)
        cell_type = "unknown"
        if isinstance(container, LoroMap):
            deep = container.get_deep_value()
            if isinstance(deep, dict):
                cell_type = str(deep.get("cellType", "unknown"))
        remaining.append({"index": i, "cell_type": cell_type})
    cell_count = len(remaining)

    if remaining:
        cell_list = ", ".join([
            f"{c['index']}: {c['cell_type']}" for c in remaining[:5]
        ])
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


@tool(
    "stop_cell",
    "Stop execution of a specific cell.",
    {
        "type": "object",
        "properties": {
            "cell_id": {"type": "string", "description": "ID of the cell to stop"},
            "cell_name": {"type": "string", "description": "Name of the cell to stop"},
            "title": {"type": "string", "description": "Title of the cell to stop"},
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

    params = {"cell_id": cell_id}

    result = await h.atomic_operation("stop_cell", params)
    if result.get("status") == "success":
        h.pending_cells.discard(cell_id)
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
async def delete_all_cells(args: dict[str, Any]) -> dict[str, Any]:  # noqa: ARG001
    try:
        notebook = await get_notebook()
        deleted_count = len(notebook.cells)
        await notebook.delete_all_cells()
    except Exception as e:
        return ok({
            "tool_name": "delete_all_cells",
            "summary": f"Failed to delete cells: {e!s}",
            "success": False,
        })

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
        "to rename the notebook or the notebook is named Untitled Layout"
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
    name = args["name"]
    print(f"[tool] rename_notebook name={name}")
    try:
        notebook = await get_notebook()
        await notebook.rename_notebook(name=name)
    except Exception as e:
        return ok({
            "tool_name": "rename_notebook",
            "summary": f"Failed to rename notebook: {e!s}",
            "success": False,
        })

    return ok({
        "tool_name": "rename_notebook",
        "success": True,
        "summary": f"Notebook renamed to '{name}'",
        "name": name,
    })


@tool(
    "create_tab",
    (
        "Create a new tab marker cell at specified position to organize cells. "
        "IMPORTANT: This inserts a new cell, shifting all subsequent cell positions down by 1."
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
    position = args["position"]
    display_name = args["display_name"]
    if position < 0:
        return ok({
            "tool_name": "create_tab",
            "summary": "Error: Position must be non-negative",
            "success": False,
        })

    print(f'[tool] create_tab pos={position} name="{display_name}"')
    try:
        notebook = await get_notebook()
        tab_id = await notebook.create_tab_marker_cell(pos=position, name=display_name)
    except Exception as e:
        return ok({
            "tool_name": "create_tab",
            "summary": f"Failed to create tab: {e!s}",
            "success": False,
        })

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
    tab_id = args["tab_id"]
    new_name = args["new_name"]

    print(f'[tool] rename_tab tab_id={tab_id} new_name="{new_name}"')
    try:
        notebook = await get_notebook()
        await notebook.rename_tab(cell_id=tab_id, new_name=new_name)
    except Exception as e:
        return ok({
            "tool_name": "rename_tab",
            "summary": f"Failed to rename tab: {e!s}",
            "success": False,
        })

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


@tool(
    "restore_checkpoint",
    (
        "Restore this notebook to one of its own previously-created checkpoints "
        "(template_version_id). Only call when the user explicitly asks to "
        "restore/revert/rollback/undo the notebook state. Confirm before restoring."
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
    template_version_id = args.get("template_version_id")
    print(f"[tool] restore_checkpoint template_version_id={template_version_id}")
    if template_version_id is None:
        return ok({
            "tool_name": "restore_checkpoint",
            "summary": "Failed to restore checkpoint: template_version_id is required",
            "success": False,
        })

    try:
        notebook = await get_notebook()
        await notebook.restore_checkpoint(template_version_id=template_version_id)
    except Exception as e:
        return ok({
            "tool_name": "restore_checkpoint",
            "summary": f"Failed to restore checkpoint: {e!s}",
            "success": False,
        })

    return ok({
        "tool_name": "restore_checkpoint",
        "success": True,
        "summary": f"Restored checkpoint to template version {template_version_id}",
        "template_version_id": template_version_id,
    })


@tool(
    "execute_code",
    (
        "Execute arbitrary Python code in the notebook kernel and return the result, "
        "stdout, stderr, and any exceptions. Use this to test imports, print values, "
        "or run simple inspection code."
    ),
    {
        "type": "object",
        "properties": {
            "code": {"type": "string", "description": "Python code to execute"}
        },
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

    print(f"[tool] execute_code: {code[:50]}...")
    result = await h.atomic_operation("execute_code", {"code": code})
    return ok({
        "tool_name": "execute_code",
        "success": True,
        "summary": "Code executed",
        "code": code,
        "stdout": result.get("stdout"),
        "stderr": result.get("stderr"),
        "exception": result.get("exception"),
    })


@tool(
    "get_global_info",
    (
        "Get rich information about a specific global variable including its type, "
        "shape, columns, dtypes, etc. Especially useful for DataFrames and AnnData objects."
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
        info = result.get("info", {})
        return ok({
            "tool_name": "get_global_info",
            "success": True,
            "summary": f"Retrieved info for global '{key}'",
            "key": key,
            "info": info,
        })

    return ok({
        "tool_name": "get_global_info",
        "success": False,
        "summary": f"Failed to get global info: {result.get('error', 'Unknown error')}",
    })


# todo(rteqs): when ClaudeCode support MCP tasks (https://modelcontextprotocol.io/seps/1686-tasks) we should make this a background task and handle tasks from SystemMessges
@tool(
    "set_widget",
    "Set a single widget value by widget key.",
    {
        "type": "object",
        "properties": {
            "key": {
                "type": "string",
                "description": (
                    "Full widget key including tf_id and widget_id in the format "
                    "<tf_id>/<widget_id>"
                ),
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
    result = await h.atomic_operation(
        "set_widget", {"key": key, "value": json.dumps(value)}
    )
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
    "get_widget",
    "Get the state of a widget by widget key (type, def, current value).",
    {
        "type": "object",
        "properties": {
            "key": {
                "type": "string",
                "description": (
                    "Full widget key including tf_id and widget_id in the format "
                    "<tf_id>/<widget_id>"
                ),
            }
        },
        "required": ["key"],
    },
)
async def get_widget(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    key = args.get("key")
    if not key:
        return ok({
            "tool_name": "get_widget",
            "summary": "Failed to get widget: Widget key is required",
            "success": False,
        })

    print(f"[tool] get_widget key={key}")
    result = await h.atomic_operation("get_widget", {"key": key})
    if result.get("status") == "success":
        state = result.get("state")
        return ok({
            "tool_name": "get_widget",
            "success": True,
            "summary": f"Retrieved widget state for: {key}",
            "key": key,
            "state": state,
        })

    return ok({
        "tool_name": "get_widget",
        "summary": f"Failed to get widget: {result.get('error', 'Unknown error')}",
        "success": False,
    })


@tool(
    "get_plot_image",
    (
        "Render a matplotlib Figure global as a webp image. Returns a base64-encoded "
        "data URL. The global must be a matplotlib Figure (or have a .figure attribute)."
    ),
    {
        "type": "object",
        "properties": {
            "key": {
                "type": "string",
                "description": "Name of the global variable holding the Figure.",
            },
            "scale": {
                "type": "number",
                "description": "DPI multiplier (base DPI 100). Default 1.0.",
            },
            "width": {
                "type": "integer",
                "description": "Output width in pixels. Default 800.",
            },
            "height": {
                "type": "integer",
                "description": "Output height in pixels. Default 600.",
            },
            "viewport": {
                "type": ["object", "null"],
                "description": (
                    "Optional axis limits: {x: [lo, hi], y: [lo, hi]}. Applied to all axes."
                ),
                "properties": {
                    "x": {"type": "array", "items": {"type": "number"}, "minItems": 2, "maxItems": 2},
                    "y": {"type": "array", "items": {"type": "number"}, "minItems": 2, "maxItems": 2},
                },
            },
        },
        "required": ["key", "scale", "width", "height"],
    },
)
async def get_plot_image(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    key = args.get("key")
    if not key:
        return ok({
            "tool_name": "get_plot_image",
            "summary": "Failed to get plot image: key is required",
            "success": False,
        })

    params = {
        "key": key,
        "scale": args.get("scale", 1.0),
        "width": args.get("width", 800),
        "height": args.get("height", 600),
        "viewport": args.get("viewport"),
    }

    print(f"[tool] get_plot_image key={key} scale={params['scale']} "
          f"{params['width']}x{params['height']} viewport={params['viewport']}")
    result = await h.atomic_operation("get_plot_image", params)
    if result.get("status") == "success":
        return ok({
            "tool_name": "get_plot_image",
            "success": True,
            "summary": f"Rendered plot image for: {key}",
            "key": key,
            "image": result.get("image"),
            "mime_type": result.get("mime_type", "image/webp"),
        })

    return ok({
        "tool_name": "get_plot_image",
        "summary": f"Failed to get plot image: {result.get('error', 'Unknown error')}",
        "success": False,
    })


@tool(
    "h5_filter_by",
    (
        "Set filters for an h5/AnnData widget. Pass the complete array of filters. "
        "Include existing filters from widget context to preserve them."
    ),
    {
        "type": "object",
        "properties": {
            "widget_key": {
                "type": "string",
                "description": (
                    "Full widget key including tf_id and widget_id in the format "
                    "<tf_id>/<widget_id>"
                ),
            },
            "label": {"type": "string", "description": "Label of the widget"},
            "filters": {
                "type": "array",
                "description": "Complete array of filters to apply",
                "items": {
                    "oneOf": [
                        {
                            "type": "object",
                            "properties": {
                                "type": {
                                    "type": "string",
                                    "enum": ["obs"],
                                    "description": "Filter by observation metadata",
                                },
                                "key": {
                                    "type": "string",
                                    "description": "The observation key to filter on",
                                },
                                "operation": {
                                    "oneOf": [
                                        {
                                            "type": "object",
                                            "properties": {
                                                "type": {
                                                    "type": "string",
                                                    "enum": ["neq"],
                                                    "description": "Not equal operation",
                                                },
                                                "value": {
                                                    "type": [
                                                        "string",
                                                        "number",
                                                        "null",
                                                    ],
                                                    "description": "Value to compare against",
                                                },
                                            },
                                            "required": ["type", "value"],
                                        },
                                        {
                                            "type": "object",
                                            "properties": {
                                                "type": {
                                                    "type": "string",
                                                    "enum": ["geq", "leq", "g", "l"],
                                                    "description": (
                                                        "Numeric comparison: geq (>=), leq (<=), "
                                                        "g (>), l (<)"
                                                    ),
                                                },
                                                "value": {
                                                    "type": "number",
                                                    "description": "Numeric value to compare against",
                                                },
                                            },
                                            "required": ["type", "value"],
                                        },
                                    ],
                                    "description": "Filter operation to apply",
                                },
                            },
                            "required": ["type", "key", "operation"],
                        },
                        {
                            "type": "object",
                            "properties": {
                                "type": {
                                    "type": "string",
                                    "enum": ["var"],
                                    "description": "Filter by variable(s) / gene(s)",
                                },
                                "keys": {
                                    "type": "array",
                                    "items": {"type": "string"},
                                    "description": "Array of variable/gene names to filter on",
                                },
                                "operation": {
                                    "type": "object",
                                    "properties": {
                                        "type": {
                                            "type": "string",
                                            "enum": ["geq", "leq", "g", "l"],
                                            "description": (
                                                "Numeric comparison: geq (>=), leq (<=), g (>), l (<)"
                                            ),
                                        },
                                        "value": {
                                            "type": "number",
                                            "description": "Numeric value to compare against",
                                        },
                                    },
                                    "required": ["type", "value"],
                                },
                            },
                            "required": ["type", "keys", "operation"],
                        },
                    ]
                },
            },
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
    result = await h.atomic_operation(
        "h5_filter_by", {"widget_key": widget_key, "filters": filters}
    )
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
    "Set an h5/AnnData widget to color by a specific observation or variable (can be multiple if for genes)",
    {
        "type": "object",
        "properties": {
            "widget_key": {
                "type": "string",
                "description": (
                    "Full widget key including tf_id and widget_id in the format "
                    "<tf_id>/<widget_id>"
                ),
            },
            "label": {"type": "string", "description": "Label of the widget"},
            "color_by": {
                "oneOf": [
                    {
                        "type": "object",
                        "properties": {
                            "type": {
                                "type": "string",
                                "enum": ["obs"],
                                "description": "Color by observation metadata",
                            },
                            "key": {
                                "type": "string",
                                "description": "The observation key to color by",
                            },
                        },
                        "required": ["type", "key"],
                    },
                    {
                        "type": "object",
                        "properties": {
                            "type": {
                                "type": "string",
                                "enum": ["var"],
                                "description": "Color by variable(s) / gene(s)",
                            },
                            "keys": {
                                "type": "array",
                                "items": {"type": "string"},
                                "description": "Array of variable/gene names to color by",
                            },
                        },
                        "required": ["type", "keys"],
                    },
                    {"type": "null"},
                ],
                "description": (
                    "Coloring configuration. Can be null to remove coloring, an obs "
                    "object to color by observation, or a var object to color by "
                    "variables like genes"
                ),
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
    result = await h.atomic_operation(
        "h5_color_by", {"widget_key": widget_key, "color_by": color_by}
    )
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
    "Refresh an h5/AnnData widget (re-render/recompute view after state changes).",
    {
        "type": "object",
        "properties": {
            "widget_key": {
                "type": "string",
                "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>",
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
    "Set the selected obsm key for an h5/AnnData widget to control which embedding is displayed.",
    {
        "type": "object",
        "properties": {
            "widget_key": {
                "type": "string",
                "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>",
            },
            "label": {"type": "string", "description": "Label of the widget"},
            "obsm_key": {
                "type": "string",
                "description": "The obsm key to use for embedding (e.g spatial, X_umap)",
            },
        },
        "required": ["widget_key", "obsm_key"],
    },
)
async def h5_set_selected_obsm_key(args: dict[str, Any]) -> dict[str, Any]:
    h = harness()
    widget_key = args.get("widget_key")
    obsm_key = args.get("obsm_key")
    print(
        f"[tool] h5_set_selected_obsm_key widget_key={widget_key} obsm_key={obsm_key}"
    )
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
    "Set a background image for an h5/AnnData widget from a user-attached file. The file must be an image type (jpg, jpeg, png, tiff).",
    {
        "type": "object",
        "properties": {
            "label": {"type": "string", "description": "Label of the widget"},
            "widget_key": {
                "type": "string",
                "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>",
            },
            "node_id": {
                "type": "string",
                "description": (
                    "The LData node ID of the image file to use as background. "
                    "This should come from files attached by the user in the chat."
                ),
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
                "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>",
            },
            "background_image_id": {
                "type": "string",
                "description": "The ID of the background image to align (typically the node_id)",
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
        f"[tool] h5_open_image_aligner widget_key={widget_key} background_image_id={background_image_id}"
    )
    params = {"widget_key": widget_key, "background_image_id": background_image_id}
    result = await h.atomic_operation("h5_open_image_aligner", params)
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
    "Reset the plotted view in an h5/AnnData widget to Plotly's autoscaled data bounds.",
    {
        "type": "object",
        "properties": {
            "label": {"type": "string", "description": "Label of the widget"},
            "widget_key": {
                "type": "string",
                "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>",
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
    "Zoom the Plotly view in an h5/AnnData widget in or out from the current camera center.",
    {
        "type": "object",
        "properties": {
            "widget_key": {
                "type": "string",
                "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>",
            },
            "label": {"type": "string", "description": "Label of the widget"},
            "direction": {
                "type": "string",
                "enum": ["in", "out"],
                "description": 'Zoom direction; use "in" to zoom closer, "out" to zoom farther',
            },
            "percentage": {
                "type": "number",
                "minimum": 0,
                "description": (
                    "Optional percentage change (e.g. 25 for ±25%); omitting uses "
                    "the default Plotly zoom factor"
                ),
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

    print(
        f"[tool] h5_zoom widget_key={widget_key} direction={direction} percentage={percentage}"
    )
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
                "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>",
            },
            "background_image_id": {
                "type": "string",
                "description": "The ID of the background image to show or hide",
            },
            "hidden": {
                "type": "boolean",
                "description": "Whether to hide (true) or show (false) the background image",
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
    "Assign selected cells to a category in a categorical observation key",
    {
        "type": "object",
        "properties": {
            "label": {"type": "string", "description": "Label of the widget"},
            "widget_key": {
                "type": "string",
                "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>",
            },
            "obs_key": {
                "type": "string",
                "description": "The existing categorical observation key to add selected cells to",
            },
            "category": {
                "type": "string",
                "description": (
                    "The category name to assign to selected cells. Will be created if it "
                    "doesn't exist in this observation key."
                ),
            },
        },
        "required": ["widget_key", "obs_key", "category"],
    },
)
async def h5_add_selected_cells_to_categorical_obs(
    args: dict[str, Any],
) -> dict[str, Any]:
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
    "Set the marker opacity for all cell markers in an h5/AnnData widget.",
    {
        "type": "object",
        "properties": {
            "label": {"type": "string", "description": "Label of the widget"},
            "widget_key": {
                "type": "string",
                "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>",
            },
            "opacity": {
                "type": "number",
                "description": "Opacity value for cell markers, between 0.1 (transparent) and 0.9 (opaque)",
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
                "description": "Full widget key including tf_id and widget_id in the format <tf_id>/<widget_id>",
            },
            "label": {"type": "string", "description": "Label of the widget"},
            "obs_key": {
                "type": "string",
                "description": "The observation column name to create or delete",
            },
            "operation": {
                "type": "string",
                "enum": ["add", "remove"],
                "description": (
                    "Whether to create a new observation column ('add') or delete an "
                    "existing one ('remove')"
                ),
            },
            "obs_type": {
                "type": "string",
                "enum": ["category", "bool", "int64", "float64"],
                "description": "Type of observation. Only for 'add' operation. Defaults to 'category'.",
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
    params = {"widget_key": widget_key, "obs_key": obs_key, "operation": operation}
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
            "package_code": {
                "type": "string",
                "description": "Multi-use package invite code.",
            },
            "package_version_id": {
                "type": "string",
                "description": "Package version ID.",
            },
            "redemption_reason": {
                "type": "string",
                "description": (
                    "Reason for redeeming the package to display in the agent history. "
                    "(e.g., `Installed X Technology Tools into Workspace`)"
                ),
            },
        },
        "required": ["package_code", "package_version_id", "redemption_reason"],
    },
)
async def redeem_package(args: dict[str, Any]) -> dict[str, Any]:
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

    package_code = str(package_code)
    package_version_id = str(package_version_id)
    print("[tool] redeem_package")
    try:
        notebook = await get_notebook()
        gql_res = await gql_query(
            query="""
                mutation RedeemPackage(
                    $argTargetAccountId: BigInt
                    $argPackageCode: String!
                    $argPackageVersionId: BigInt!
                    $argAllowSharing: Boolean
                ) {
                    redeemPackage(
                        input: {
                            argTargetAccountId: $argTargetAccountId
                            argPackageCode: $argPackageCode
                            argPackageVersionId: $argPackageVersionId
                            argAllowSharing: $argAllowSharing
                        }
                    ) {
                        clientMutationId
                        result {
                            resPackageId
                            resPackageVersionId
                            resRedeemedInAccountId
                        }
                    }
                }
            """,
            variables={
                "argTargetAccountId": notebook.owner_id,
                "argPackageCode": package_code,
                "argPackageVersionId": package_version_id,
                "argAllowSharing": True,
            },
            auth=auth_token_sdk,
        )
        validate(gql_res, RedeemPackageRes)
    except Exception as e:
        return ok({
            "tool_name": "redeem_package",
            "success": False,
            "summary": f"Failed to redeem package: {e!s}",
        })

    return ok({
        "tool_name": "redeem_package",
        "success": True,
        "summary": f"Redeemed package {package_code} (version {package_version_id})",
        "package_code": package_code,
        "package_version_id": package_version_id,
        "redemption_reason": redemption_reason,
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
                "description": (
                    "Optional full widget key including tf_id and widget_id in the format "
                    "<tf_id>/<widget_id> for keywords related to a specific widget"
                ),
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
    params = {"keyword": keyword}
    if widget_key is not None:
        params["widget_key"] = widget_key

    result = await h.atomic_operation("smart_ui_spotlight", params)
    if result.get("status") == "success":
        res = {
            "tool_name": "smart_ui_spotlight",
            "success": True,
            "summary": f"Highlighted UI element: {keyword}",
            "keyword": keyword,
        }
        if widget_key is not None:
            res["widget_key"] = widget_key
        if widget_label is not None:
            res["widget_label"] = widget_label
        return ok(res)

    return ok({
        "tool_name": "smart_ui_spotlight",
        "success": False,
        "summary": f"Failed to highlight UI element: {result.get('error', 'Unknown error')}",
    })


@tool(
    "update_plan",
    "Update the agent's plan.",
    {
        "type": "object",
        "properties": {
            "plan": {
                "type": "array",
                "description": "List of plan items. This should be the complete current plan state.",
                "items": {
                    "type": "object",
                    "properties": {
                        "id": {
                            "type": "string",
                            "description": "Unique step identifier",
                        },
                        "description": {
                            "type": "string",
                            "description": "What this step does",
                        },
                        "status": {
                            "type": "string",
                            "enum": ["todo", "in_progress", "done", "cancelled"],
                            "description": "Current status",
                        },
                    },
                    "required": ["id", "description", "status"],
                },
            },
            "plan_diff": {
                "type": "array",
                "description": "List of plan diff items - only the steps that changed in this update.",
                "items": {
                    "type": "object",
                    "properties": {
                        "id": {
                            "type": "string",
                            "description": "Unique step identifier",
                        },
                        "action": {
                            "type": "string",
                            "enum": ["add", "update", "complete", "remove"],
                            "description": "What action was taken on this step",
                        },
                    },
                    "required": ["id", "action"],
                },
            },
            "plan_update_overview": {
                "type": "string",
                "description": (
                    "Short title overview of what changed. E.g. 'Added QC steps' or "
                    "'Completed step 2, step 3 now in progress'"
                ),
            },
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

        plan = {"steps": plan_items}
        h.current_plan = plan

        print(f"[tool] update_plan: {plan_update_overview}")
        for item in plan_items:
            print(
                f"  [{item.get('status')}] {item.get('id')}: {item.get('description')}"
            )
        for diff in plan_diff:
            print(f"  diff: [{diff.get('action')}] {diff.get('id')}")

        return ok({
            "tool_name": "update_plan",
            "success": True,
            "summary": plan_update_overview or "Plan updated",
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
    "Submit the user-facing response with next_status, questions, and summary. Call this at the end of every loop.",
    {
        "type": "object",
        "properties": {
            "summary": {
                "type": "string",
                "description": "User-facing progress, responses, or next step. Use markdown bullets if needed.",
            },
            "next_status": {
                "type": "string",
                "description": "What the agent will do next",
                "enum": [
                    "executing",
                    "fixing",
                    "thinking",
                    "awaiting_user_response",
                    "awaiting_cell_execution",
                    "awaiting_user_widget_input",
                    "done",
                ],
            },
            "expected_widgets": {
                "type": "array",
                "items": {"type": "string"},
                "description": (
                    "Optional list of full widget keys (<tf_id>/<widget_id>) to await "
                    "when next_status is 'awaiting_user_widget_input'"
                ),
            },
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

        next_status = args.get("next_status", "done")
        if not isinstance(next_status, str) or next_status not in {
            "executing",
            "fixing",
            "thinking",
            "awaiting_user_response",
            "awaiting_cell_execution",
            "awaiting_user_widget_input",
            "done",
        }:
            print(f"[agent] Invalid next_status: {next_status}")
            return ok({
                "tool_name": "submit_response",
                "success": False,
                "message": "Please provide a valid next_status",
            })

        expected_widgets = args.get("expected_widgets", [])

        print("[tool] submit_response called with:")
        print(f"  - next_status: {next_status}")
        print(f"  - summary: {summary}")
        print(f"  - expected_widgets: {expected_widgets}")

        h.current_status = next_status
        if next_status == "awaiting_user_widget_input":
            h.pending_widgets = {str(k): None for k in expected_widgets}

        await h.set_agent_status(next_status)
        return ok({
            "tool_name": "submit_response",
            "success": True,
            "summary": "Response submitted successfully",
        })
    except Exception as e:
        print(f"[tool] submit_response error: {e!s}")
        import traceback

        traceback.print_exc()
        return ok({
            "tool_name": "submit_response",
            "success": False,
            "summary": f"Error submitting response: {e!s}",
        })


all_tools = [
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
    update_plan,
    set_widget,
    get_widget,
    get_plot_image,
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
    redeem_package,
    smart_ui_spotlight,
    submit_response,
]

agent_tools_mcp = create_sdk_mcp_server(name=MCP_SERVER_NAME, tools=all_tools)


async def can_use_tool(
    tool_name: str, input_data: dict[str, Any], context: ToolPermissionContext
) -> PermissionResultAllow | PermissionResultDeny:
    if tool_name == "AskUserQuestion":
        h = harness()

        await h.pending_question_event.wait()
        tx_id = h.pending_question_tx_id
        assert tx_id is not None

        await h.set_agent_status("awaiting_user_response")

        try:
            result = await h.wait_for_operation(tx_id, timeout=300)
        except Exception as e:
            print(f"[agent] AskUserQuestion failed: {e!s}")
            return PermissionResultDeny(
                message=f"Failed to get user response: {e!s}", interrupt=True
            )
        finally:
            h.pending_question_event.clear()

        if result.get("status") == "success":
            if result.get("skip") is True:
                return PermissionResultDeny(
                    message="User skipped question", interrupt=True
                )

            questions = input_data.get("questions")
            answers = result.get("answers", {})

            qa_content = {
                "type": "answers",
                "content": [{"question": q, "answer": a} for q, a in answers.items()],
            }
            await h._insert_history(payload={"content": qa_content})

            return PermissionResultAllow(
                updated_input={"questions": questions, "answers": answers}
            )

        return PermissionResultDeny(
            message=f"Failed to get user response: {result.get('error', 'Unknown error')}",
            interrupt=True,
        )

    return PermissionResultAllow(updated_input=input_data)
