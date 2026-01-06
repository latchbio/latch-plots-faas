"""Endpoints for headless browser management."""

import orjson
from latch_asgi.context.http import Context, HandlerResult

from ..entrypoint import (
    headless_browser,
    plots_ctx_manager,
    start_headless_browser_if_enabled,
)


async def spawn_headless_browser(ctx: Context) -> HandlerResult:
    """
    Trigger headless browser spawn.
    Called after frontend sets headless_mode=true via GraphQL.
    
    Returns:
        Status of the spawn operation.
    """
    global headless_browser

    if headless_browser is not None:
        return orjson.dumps({"status": "already_running"})

    notebook_id = plots_ctx_manager.notebook_id
    if notebook_id is None:
        return orjson.dumps({"status": "error", "message": "No notebook ID available"})

    try:
        await start_headless_browser_if_enabled(str(notebook_id))
        
        if headless_browser is not None:
            return orjson.dumps({"status": "spawned"})
        else:
            return orjson.dumps({"status": "not_enabled", "message": "headless_mode not enabled in database"})
    except Exception as e:
        return orjson.dumps({"status": "error", "message": str(e)})


async def stop_headless_browser_endpoint(ctx: Context) -> HandlerResult:
    """
    Stop the headless browser if running.
    
    Returns:
        Status of the stop operation.
    """
    from ..entrypoint import stop_headless_browser

    if headless_browser is None:
        return orjson.dumps({"status": "not_running"})

    try:
        await stop_headless_browser()
        return orjson.dumps({"status": "stopped"})
    except Exception as e:
        return orjson.dumps({"status": "error", "message": str(e)})
