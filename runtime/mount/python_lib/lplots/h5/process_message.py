import time
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
            if name is not None and name.endswith(".parquet"):
                transcript_path = f
                break

        if transcript_path is None:
            return {
                "type": "h5",
                "data_type": "transcripts",
                "key": widget_session_key,
                "value": {
                    "error": "Invalid message -- no parquet files found in spatial directory",
                },
            }

        duckdb_table_name = f"h5spatial_{transcript_path.version_id()}"

        table_exists = _inject.kernel.duckdb.execute("""
            SELECT EXISTS (
                SELECT 1
                FROM information_schema.tables
                WHERE table_name = ?
            )
        """, [duckdb_table_name]).fetchone()
        assert table_exists is not None
        table_exists = table_exists[0]

        create_table_time = 0
        if not table_exists:
            start_time = time.time()
            _inject.kernel.duckdb.execute(f"""
                CREATE TABLE {duckdb_table_name} (
                    fov INTEGER,
                    cell_id INTEGER,
                    target TEXT,
                    cell_comp TEXT,
                    global_x DOUBLE,
                    global_y DOUBLE
                )
            """)
            presigned_url = await _inject.kernel.get_presigned_url(transcript_path.path)
            _inject.kernel.duckdb.execute(f"""
                COPY {duckdb_table_name}
                FROM '{presigned_url}'
                (FORMAT 'PARQUET')
            """)
            create_table_time = round(time.time() - start_time, 2)

        return await process_spatial_request(msg, widget_session_key, duckdb_table_name, create_table_time)

    raise ValueError(f"Invalid H5 viewer message: {msg}")
