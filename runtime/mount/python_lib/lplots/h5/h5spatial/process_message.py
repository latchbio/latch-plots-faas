import time
from typing import Any

import duckdb
from duckdb import (
    ColumnExpression,
    ConstantExpression,
)

from ... import _inject


def get_spatial_sample(
    conn: duckdb.DuckDBPyConnection,
    table_name: str,
    x_min: float,
    y_min: float,
    x_max: float,
    y_max: float,
    max_transcripts: int = 100000,
    genes_to_fetch: list[str] | None = None,
) -> tuple[duckdb.DuckDBPyRelation, int, int]:
    table_rel = conn.sql(f"select * from {table_name}")  # noqa: S608

    spatial_filter = (
        (ColumnExpression("global_x") >= ConstantExpression(x_min)) &
        (ColumnExpression("global_x") <= ConstantExpression(x_max)) &
        (ColumnExpression("global_y") >= ConstantExpression(y_min)) &
        (ColumnExpression("global_y") <= ConstantExpression(y_max))
    )

    if genes_to_fetch is not None and len(genes_to_fetch) > 0:
        gene_expressions = [ConstantExpression(gene) for gene in genes_to_fetch]
        gene_filter = ColumnExpression("target").isin(*gene_expressions)
        filter_condition = spatial_filter & gene_filter
    else:
        filter_condition = spatial_filter

    points_in_scope_rel = table_rel.filter(filter_condition).aggregate(
        "count(*) as cnt"
    )
    points_in_scope_result = points_in_scope_rel.fetchone()
    points_in_scope = 0 if points_in_scope_result is None else points_in_scope_result[0]

    total_points_rel = table_rel.aggregate(
        "count(*)"
    )
    total_points_result = total_points_rel.fetchone()
    total_points = 0 if total_points_result is None else total_points_result[0]

    sampled_rel = (
        table_rel
        .filter(filter_condition)
        .order("random()")
        .limit(max_transcripts)
    )

    return sampled_rel, points_in_scope, total_points


async def process_spatial_request(  # noqa: RUF029
    msg: dict[str, Any],
    widget_session_key: str,
    duckdb_table_name: str,
    create_table_time: float,
) -> dict[str, Any]:
    if "op" not in msg or msg["op"] not in {"get"}:
        return {
            "type": "h5",
            "key": widget_session_key,
            "data_type": "transcripts",
            "value": {
                "error": f"Invalid transcripts operation: {msg.get('op', '`op` key missing from message')}",
            },
        }

    op = msg["op"]

    max_transcripts = int(msg.get("max_transcripts", 100000))

    genes_to_fetch = msg.get("genes")

    if op == "get":
        start_time = time.time()
        sampled_data, points_in_scope, total_points = get_spatial_sample(
            _inject.kernel.duckdb,
            duckdb_table_name,
            x_min=float(msg["x_min"]),
            y_min=float(msg["y_min"]),
            x_max=float(msg["x_max"]),
            y_max=float(msg["y_max"]),
            max_transcripts=max_transcripts,
            genes_to_fetch=genes_to_fetch,
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
                    "fetched_for_max_transcripts": max_transcripts,
                    "fetched_for_x_min": float(msg["x_min"]),
                    "fetched_for_x_max": float(msg["x_max"]),
                    "fetched_for_y_min": float(msg["y_min"]),
                    "fetched_for_y_max": float(msg["y_max"]),
                    "fetched_for_genes": genes_to_fetch,
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


def get_boundary_sample(
    conn: duckdb.DuckDBPyConnection,
    table_name: str,
    x_min: float, y_min: float, x_max: float, y_max: float,
    max_boundaries: int = 100_000,
) -> duckdb.DuckDBPyRelation:
    return conn.sql(f"""
        select
            ID,
            EntityID,
            case
                when st_geometrytype(Geometry) = 'POLYGON' then 
                    st_astext(st_exteriorring(Geometry))
                when st_geometrytype(Geometry) = 'LINESTRING' then
                    st_astext(Geometry)
                else
                    st_astext(Geometry)
            end as coords
        from
            {table_name}
        where
            ST_Intersects(Geometry, ST_GeomFromText('POLYGON(({x_min} {y_min}, {x_max} {y_min}, {x_max} {y_max}, {x_min} {y_max}, {x_min} {y_min}))'))
        order by random()
        limit {max_boundaries}
    """)  # noqa: S608


async def process_boundaries_request(  # noqa: RUF029
    msg: dict[str, Any],
    widget_session_key: str,
    duckdb_table_name: str,
    create_table_time: float,
) -> dict[str, Any]:
    if "op" not in msg or msg["op"] not in {"get"}:
        return {
            "type": "h5",
            "key": widget_session_key,
            "data_type": "boundaries",
            "value": {
                "error": f"Invalid boundaries operation: {msg.get('op', '`op` key missing from message')}",
            },
        }

    op = msg["op"]

    max_boundaries = int(msg.get("max_boundaries", 100_000))

    if op == "get":
        start_time = time.time()
        sampled_data = get_boundary_sample(
            _inject.kernel.duckdb,
            duckdb_table_name,
            x_min=float(msg["x_min"]),
            y_min=float(msg["y_min"]),
            x_max=float(msg["x_max"]),
            y_max=float(msg["y_max"]),
            max_boundaries=max_boundaries,
        )

        columns = sampled_data.columns
        data = sampled_data.fetchall()

        entity_ids = []
        ids = []
        boundaries = []

        for row in data:
            row_id: str = row[0]
            entity_id: str = row[1]
            wkt_coords: str | None = row[2]

            entity_ids.append(entity_id)
            ids.append(row_id)

            if wkt_coords is None:
                boundaries.append(None)
                continue

            if wkt_coords.startswith("LINESTRING("):
                coords_str = wkt_coords[11:-1]
            elif wkt_coords.startswith("POLYGON(("):
                coords_str = wkt_coords[9:-2]  # Remove "POLYGON((" and "))"
            else:
                coords_str = wkt_coords

            coord_pairs = []
            for pair in coords_str.split(","):
                x, y = pair.strip().split()
                coord_pairs.append([float(x), float(y)])

            boundaries.append(coord_pairs)

        return {
            "type": "h5",
            "op": op,
            "data_type": "boundaries",
            "key": widget_session_key,
            "value": {
                "data": {
                    "columns": columns,
                    "entity_ids": entity_ids,
                    "ids": ids,
                    "boundaries": boundaries,
                    "time_taken": round(time.time() - start_time, 2),
                    "create_table_time": create_table_time,
                    "fetched_for_max_boundaries": max_boundaries,
                    "fetched_for_x_min": float(msg["x_min"]),
                    "fetched_for_x_max": float(msg["x_max"]),
                    "fetched_for_y_min": float(msg["y_min"]),
                    "fetched_for_y_max": float(msg["y_max"]),
                }
            },
        }

    return {
        "type": "h5",
        "data_type": "boundaries",
        "key": widget_session_key,
        "value": {
            "error": f"Unsupported boundaries operation: {op}",
        },
    }
