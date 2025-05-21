from typing import Any

from lplots.widgets.h5 import H5AD, H5Spatial
from runtime.mount.python_lib.lplots.h5.h5ad.process_message import process_h5ad_request
from runtime.mount.python_lib.lplots.h5.h5spatial.process_message import (
    process_h5spatial_request,
)
from runtime.mount.python_lib.lplots.h5.utils import auto_install
from runtime.mount.utils import get_presigned_url

from .. import _inject

ad = auto_install.ad


async def handle_h5_widget_message(
    msg: dict[str, Any],
) -> dict[str, Any]:
    if msg["type"] != "ann_data" or "key" not in msg or "state" not in msg:
        return {
            "type": "ann_data",
            "key": None,
            "value": {
                "error": "Invalid message -- missing `key` or `state`",
            },
        }

    widget_session_key = msg["key"]
    widget_state: H5AD | H5Spatial | Any = msg["state"]["data"]

    if isinstance(widget_state, H5AD):
        obj_id = widget_state.obj_id
        if obj_id is None or obj_id not in _inject.kernel.ann_data_objects:
            return {
                "type": "ann_data",
                "key": widget_session_key,
                "value": {
                    "error": f"AnnData object with ID {obj_id} not found in cache",
                },
            }

        adata: ad.AnnData = _inject.kernel.ann_data_objects[obj_id]

        return await process_h5ad_request(msg, widget_session_key, adata, obj_id)

    if isinstance(widget_state, H5Spatial):
        duckdb_table_name = f"h5_spatial_{widget_state.transcript_path.version_id()}"

        table_exists = _inject.kernel.duckdb.execute("""
            SELECT EXISTS (
                SELECT 1
                FROM information_schema.tables
                WHERE table_name = ?
            )
        """, [duckdb_table_name]).fetchone()[0]

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
            presigned_url = await get_presigned_url(widget_state.transcript_path.path)
            _inject.kernel.duckdb.execute(f"""
                COPY {duckdb_table_name}
                FROM '{presigned_url}'
                (FORMAT 'PARQUET')
            """)

        return await process_h5spatial_request(msg, widget_session_key, duckdb_table_name)

    raise ValueError(f"Invalid H5 viewer message: {msg}")
