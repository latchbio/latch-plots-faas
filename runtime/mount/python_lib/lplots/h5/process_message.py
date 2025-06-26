import time
from collections.abc import Awaitable, Callable
from pathlib import Path
from typing import Any, Literal

from duckdb import (
    ColumnExpression,
    ConstantExpression,
)
from latch.ldata.path import LPath

from lplots.h5.h5ad.process_message import process_h5ad_request
from lplots.h5.h5spatial.process_message import (
    process_boundaries_request,
    process_spatial_request,
)
from lplots.h5.utils import auto_install

from .. import _inject

ad = auto_install.ad


async def handle_h5_widget_message(
    msg: dict[str, Any],
    send: Callable[[object], Awaitable[None]]
) -> dict[str, Any]:
    if msg["type"] != "h5" or "key" not in msg or "state" not in msg:
        return {
            "type": "h5",
            "key": None,
            "value": {
                "error": "Invalid message -- missing `key` or `state`",
            },
        }

    widget_session_key: str = msg["key"]
    data_type: Literal["h5ad", "transcripts", "boundaries"] | Any = msg.get("data_type", "h5ad")
    widget_state: dict[str, Any] = msg["state"]

    if data_type == "h5ad":
        obj_id = widget_state["obj_id"]
        if obj_id is None or obj_id not in _inject.kernel.ann_data_objects:
            return {
                "type": "h5",
                "data_type": "h5ad",
                "key": widget_session_key,
                "value": {
                    "error": f"AnnData object with ID {obj_id} not found in cache",
                },
            }

        adata: ad.AnnData = _inject.kernel.ann_data_objects[obj_id]

        return await process_h5ad_request(msg, widget_session_key, adata,
                                          obj_id, send)

    if data_type in {"transcripts", "boundaries"}:
        if "spatial_dir" not in widget_state:
            return {
                "type": "h5",
                "data_type": data_type,
                "key": widget_session_key,
                "value": {
                    "error": "Invalid message -- missing `spatial_dir`",
                },
            }

        spatial_dir = LPath(widget_state["spatial_dir"]["path"])

        duckdb_file_path = None
        for f in spatial_dir.iterdir():
            name = f.name()
            if name == "transcripts.duckdb":
                duckdb_file_path = f
                break

        if duckdb_file_path is None:
            return {
                "type": "h5",
                "data_type": data_type,
                "key": widget_session_key,
                "value": {
                    "error": f"Invalid message -- no duckdb files found in spatial directory for {data_type}",
                },
            }

        sanitized_widget_session_key = widget_session_key.replace(":", "_").replace("/", "_").replace("-", "_")
        schema_name = f"{data_type}_{sanitized_widget_session_key}"
        table_name = "final_transcripts" if data_type == "transcripts" else "cell_boundaries"

        db_databases_rel = _inject.kernel.duckdb.sql("select * from duckdb_databases()")
        db_attached_rel = db_databases_rel.filter(
            ColumnExpression("database_name") == ConstantExpression(schema_name)
        ).aggregate("count(*) > 0 as is_attached")

        db_attached_result = db_attached_rel.fetchone()
        db_attached = db_attached_result[0] if db_attached_result else False

        create_table_time = 0
        if not db_attached:
            start_time = time.time()

            local_duckdb_file_path = Path(f"/tmp/{data_type}_{sanitized_widget_session_key}.duckdb")  # noqa: S108
            duckdb_file_path.download(local_duckdb_file_path, cache=True)

            _inject.kernel.duckdb.execute(f"attach '{local_duckdb_file_path}' as {schema_name} (read_only)")
            create_table_time = round(time.time() - start_time, 2)

        duckdb_table_name = f"{schema_name}.{table_name}"

        if data_type == "transcripts":
            return await process_spatial_request(msg, widget_session_key, duckdb_table_name, create_table_time)

        if data_type == "boundaries":
            return await process_boundaries_request(msg, widget_session_key, duckdb_table_name, create_table_time)

    raise ValueError(f"Invalid H5 viewer message: {msg}")
