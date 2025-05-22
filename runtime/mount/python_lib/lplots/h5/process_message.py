from typing import Any

from latch.ldata.path import LPath

from lplots.h5.h5ad.process_message import process_h5ad_request
from lplots.h5.h5spatial.process_message import (
    process_h5spatial_request,
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
    widget_state: dict[str, Any] = msg["state"]["data"]

    if "obj_id" in widget_state:
        obj_id = widget_state["obj_id"]
        if obj_id is None or obj_id not in _inject.kernel.ann_data_objects:
            return {
                "type": "h5",
                "key": widget_session_key,
                "value": {
                    "error": f"AnnData object with ID {obj_id} not found in cache",
                },
            }

        adata: ad.AnnData = _inject.kernel.ann_data_objects[obj_id]

        return await process_h5ad_request(msg, widget_session_key, adata, obj_id)

    if "transcript_path" in widget_state:
        transcript_path = LPath(widget_state["transcript_path"]["path"])
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

        if not table_exists:
            _inject.kernel.duckdb.execute(f"""
                CREATE TABLE {duckdb_table_name} (
                    fov INTEGER,
                    cell_id INTEGER,
                    x DOUBLE,
                    y DOUBLE,
                    target TEXT,
                    cell_comp TEXT
                )
            """)
            presigned_url = await _inject.kernel.get_presigned_url(transcript_path.path)
            _inject.kernel.duckdb.execute(f"""
                COPY {duckdb_table_name}
                FROM '{presigned_url}'
                (FORMAT 'PARQUET')
            """)

        return await process_h5spatial_request(msg, widget_session_key, duckdb_table_name)

    raise ValueError(f"Invalid H5 viewer message: {msg}")
