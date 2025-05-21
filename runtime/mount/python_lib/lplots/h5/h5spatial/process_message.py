from typing import Any

import duckdb
from duckdb import DuckDBPyRelation

from ... import _inject


def get_spatial_sample(
    conn: duckdb.DuckDBPyConnection,
    table_name: str,
    x_min: float,
    y_min: float,
    x_max: float,
    y_max: float,
    max_points: int = 100000,
) -> DuckDBPyRelation:
    count_result = conn.sql(f"""
        SELECT COUNT(*) as cnt
        FROM {table_name}
        WHERE x >= {x_min} AND x <= {x_max}
            AND y >= {y_min} AND y <= {y_max}
    """).fetchone()  # noqa: S608

    total_points = 0 if count_result is None else count_result[0]

    if total_points <= max_points:
        return conn.sql(f"""
            SELECT *
            FROM {table_name}
            WHERE x >= {x_min} AND x <= {x_max}
                AND y >= {y_min} AND y <= {y_max}
        """)  # noqa: S608

    max_occupancy = 2
    total_cells = max_points // max_occupancy

    viewport_width = x_max - x_min
    viewport_height = y_max - y_min

    cells_x = int((total_cells * viewport_width / viewport_height) ** 0.5)
    cells_y = total_cells // cells_x

    spatial_data = conn.sql(f"""
        SELECT
            *,
            row_number() OVER (
                PARTITION BY
                    floor(({cells_x} * (x - {x_min}) / {viewport_width}) +
                    {cells_x} * floor({cells_y} * (y - {y_min}) / {viewport_height}))
            ) as row_num
        FROM {table_name}
        WHERE x >= {x_min} AND x <= {x_max}
            AND y >= {y_min} AND y <= {y_max}
    """)  # noqa: S608

    filtered_data = spatial_data.filter(f"row_num <= {max_occupancy}")
    return filtered_data.order("random()").limit(max_points)


async def process_h5spatial_request(  # noqa: RUF029
    msg: dict[str, Any],
    widget_session_key: str,
    duckdb_table_name: str,
) -> dict[str, Any]:
    if "op" not in msg or msg["op"] not in {"init_data"}:  # noqa: FURB171
        return {
            "type": "ann_data",
            "key": widget_session_key,
            "value": {
                "error": f"Invalid h5spatial operation: {msg.get('op', '`op` key missing from message')}",
            },
        }

    op = msg["op"]

    if op == "init_data":
        sampled_data = get_spatial_sample(
            _inject.kernel.duckdb,
            duckdb_table_name,
            x_min=float(msg["x_min"]),
            y_min=float(msg["y_min"]),
            x_max=float(msg["x_max"]),
            y_max=float(msg["y_max"]),
        )

        columns = sampled_data.columns
        data = sampled_data.fetchall()

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_session_key,
            "value": {
                "data": {
                    "columns": columns,
                    "transcripts": data,
                }
            },
        }

    return {
        "type": "ann_data",
        "key": widget_session_key,
        "value": {
            "error": f"Unsupported operation: {op}",
        },
    }
