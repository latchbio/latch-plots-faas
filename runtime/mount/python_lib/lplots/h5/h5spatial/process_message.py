import time
from typing import Any

import duckdb

from ... import _inject


def get_spatial_sample(
    conn: duckdb.DuckDBPyConnection,
    table_name: str,
    x_min: float,
    y_min: float,
    x_max: float,
    y_max: float,
    max_points: int = 100000,
) -> tuple[duckdb.DuckDBPyRelation, int, int]:
    points_in_scope_q = conn.sql(f"""
        SELECT COUNT(*) AS cnt
        FROM {table_name}
        WHERE global_x >= {x_min} AND global_x <= {x_max}
          AND global_y >= {y_min} AND global_y <= {y_max}
    """).fetchone()  # noqa: S608
    points_in_scope = 0 if points_in_scope_q is None else points_in_scope_q[0]

    total_points_q = conn.sql(f"""SELECT COUNT(*) AS cnt FROM {table_name}""").fetchone()  # noqa: S608
    total_points = 0 if total_points_q is None else total_points_q[0]

    sampled_rel = conn.sql(f"""
        SELECT *
        FROM (
            SELECT *
            FROM {table_name}
            WHERE global_x BETWEEN {x_min} AND {x_max}
              AND global_y BETWEEN {y_min} AND {y_max}
        ) USING SAMPLE reservoir({max_points} ROWS)
    """)  # noqa: S608

    return sampled_rel, points_in_scope, total_points


async def process_spatial_request(  # noqa: RUF029
    msg: dict[str, Any],
    widget_session_key: str,
    duckdb_table_name: str,
    create_table_time: float,
) -> dict[str, Any]:
    if "op" not in msg or msg["op"] not in {"get"}:  # noqa: FURB171
        return {
            "type": "h5",
            "key": widget_session_key,
            "data_type": "transcripts",
            "value": {
                "error": f"Invalid transcripts operation: {msg.get('op', '`op` key missing from message')}",
            },
        }

    op = msg["op"]

    max_points = int(msg.get("max_points", 100000))

    if op == "get":
        start_time = time.time()
        sampled_data, points_in_scope, total_points = get_spatial_sample(
            _inject.kernel.duckdb,
            duckdb_table_name,
            x_min=float(msg["x_min"]),
            y_min=float(msg["y_min"]),
            x_max=float(msg["x_max"]),
            y_max=float(msg["y_max"]),
            max_points=max_points,
        )

        columns = sampled_data.columns
        data = sampled_data.fetchall()

        return {
            "type": "h5",
            "op": op,
            "data_type": "transcripts",
            "key": widget_session_key,
            "value": {
                "data": {
                    "columns": columns,
                    "transcripts": data,
                    "points_in_scope": points_in_scope,
                    "total_points": total_points,
                    "time_taken": round(time.time() - start_time, 2),
                    "create_table_time": create_table_time,
                    "fetched_for_max_points": max_points,
                    "fetched_for_x_min": float(msg["x_min"]),
                    "fetched_for_x_max": float(msg["x_max"]),
                    "fetched_for_y_min": float(msg["y_min"]),
                    "fetched_for_y_max": float(msg["y_max"]),
                }
            },
        }

    return {
        "type": "h5",
        "data_type": "transcripts",
        "key": widget_session_key,
        "value": {
            "error": f"Unsupported transcripts operation: {op}",
        },
    }
