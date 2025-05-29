import time
from pathlib import Path
from typing import Any, Literal

from latch.ldata.path import LPath

from lplots.h5.h5ad.process_message import process_h5ad_request
from lplots.h5.h5spatial.process_message import (
    process_spatial_request,
)
from lplots.h5.utils import auto_install

from .. import _inject

ad = auto_install.ad


async def handle_h5_widget_message(
    msg: dict[str, Any],
) -> dict[str, Any]:
    if msg["type"] != "h5" or "key" not in msg or "state" not in msg:
        return {
            "type": "h5",
            "key": None,
            "value": {
                "error": "Invalid message -- missing `key` or `state`",
            },
        }

    widget_session_key = msg["key"]
    data_type: Literal["h5ad", "transcripts"] | Any = msg.get("data_type", "h5ad")
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

        return await process_h5ad_request(msg, widget_session_key, adata, obj_id)

    if data_type == "transcripts":
        if "spatial_dir" not in widget_state:
            return {
                "type": "h5",
                "data_type": "transcripts",
                "key": widget_session_key,
                "value": {
                    "error": "Invalid message -- missing `spatial_dir`",
                },
            }

        spatial_dir = LPath(widget_state["spatial_dir"]["path"])
        transcript_path = None
        for f in spatial_dir.iterdir():
            name = f.name()
            if name is not None and name.endswith(".duckdb"):
                transcript_path = f
                break

        if transcript_path is None:
            return {
                "type": "h5",
                "data_type": "transcripts",
                "key": widget_session_key,
                "value": {
                    "error": "Invalid message -- no duckdb files found in spatial directory",
                },
            }

        schema_name = "transcripts"
        table_name = "final_transcripts"

        db_attached = _inject.kernel.duckdb.execute("""
            select
                count(*) > 0
            from
                duckdb_databases()
            where
                database_name = ?
        """, [schema_name]).fetchone()
        db_attached = db_attached[0]

        create_table_time = 0
        if not db_attached:
            start_time = time.time()
            local_transcript_path = Path("/tmp/transcripts.duckdb")  # noqa: S108
            transcript_path.download(local_transcript_path, cache=True)

            _inject.kernel.duckdb.execute(f"attach '{local_transcript_path}' as transcripts (read_only)")
            create_table_time = round(time.time() - start_time, 2)

        duckdb_table_name = f"{schema_name}.{table_name}"
        return await process_spatial_request(msg, widget_session_key, duckdb_table_name, create_table_time)

    raise ValueError(f"Invalid H5 viewer message: {msg}")
