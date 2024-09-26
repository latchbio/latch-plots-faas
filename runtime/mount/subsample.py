import datetime
import sys
import traceback
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, NotRequired, TypedDict

from duckdb import (
    ColumnExpression,
    ConstantExpression,
    DuckDBPyConnection,
    DuckDBPyRelation,
)
from duckdb import connect as duckdb_connect
from pandas import DataFrame
from plotly.basedatatypes import BaseFigure

sys.path.append(str(Path(__file__).parent.absolute()))
from utils import auth_token_sdk, get_presigned_url, gql_query


def initialize_duckdb() -> DuckDBPyConnection:
    conn = duckdb_connect(database=":memory:plots-faas", read_only=False)

    conn.execute(
        """
        create table plots_faas_catalog (
            name string not null unique,
            generation integer,
            last_modified_time timestamp
        )
        """
    )
    return conn


class Trace(TypedDict):
    x: str
    y: str
    color_by: NotRequired[str]
    error_bar: NotRequired[str]
    marker_size: NotRequired[str]


class PlotConfig(TypedDict):
    traces: list[Trace]
    custom_data: NotRequired[list[str]]

    facet: NotRequired[str]

    xrange: NotRequired[tuple[float | int, float | int]]
    yrange: NotRequired[tuple[float | int, float | int]]
    height_px: NotRequired[float]
    width_px: NotRequired[float]


async def check_generation(
    conn: DuckDBPyConnection,
    table_name: str,
    source_type: Literal["kernel", "ldata"] = "kernel",
    cur_gen: int | None = None,
) -> tuple[bool, datetime.datetime | None]:
    last_modified_time = None
    if source_type == "ldata":
        try:
            resp = await gql_query(
                auth=auth_token_sdk,
                query="""
                    query GetLDataLastModifiedTime($ldataNodeId: BigInt!) {
                        ldataNodeEvents(
                            filter: {
                            type: { equalTo: INGRESS }
                            ldataNodeId: { equalTo: $ldataNodeId }
                            }
                            orderBy: TIME_DESC
                            first: 1
                        ) {
                            nodes {
                              id
                              time
                            }
                        }
                    }
                """,
                variables={"ldataNodeId": int(table_name.split("_")[1])},
            )
            data = resp["data"]["ldataNodeEvents"]["nodes"]
            if len(data) > 0:
                # todo(rteqs): check if this is right
                last_modified_time = datetime.datetime.fromisoformat(data[0]["time"])
        except Exception:
            traceback.print_exc()
            raise

    last_modified_time = datetime.datetime.now(tz=datetime.UTC)

    duckdb_gen = (
        conn.table("plots_faas_catalog")
        .filter(ColumnExpression("name") == ConstantExpression(table_name))
        .project(
            (ColumnExpression("generation") == ConstantExpression(cur_gen))
            or (
                ColumnExpression("last_modified_time")
                > ConstantExpression(last_modified_time)
            )
        )
        .fetchone()
    )

    if duckdb_gen is None:
        return False, None

    if not isinstance(duckdb_gen[0], bool):
        raise TypeError(f"Unexpected generation value: {duckdb_gen}")

    if TYPE_CHECKING:
        assert isinstance(duckdb_gen[0], bool)

    return duckdb_gen[0], last_modified_time


# todo(rteqs): marker size
def downsample(
    conn: DuckDBPyConnection, table_name: str, config: PlotConfig
) -> list[DuckDBPyRelation]:
    custom_data = config.get("custom_data")
    custom_data_str = (
        ", ".join(custom_data)
        if custom_data is not None and len(custom_data) > 0
        else None
    )

    facet = config.get("facet")
    marker_size_axis = config.get("marker_size_axis")
    width_px = config.get("width_px")
    height_px = config.get("height_px")

    facet_and_color_by: set[str] = set()
    if facet is not None:
        facet_and_color_by.add(facet)
    for trace in config["traces"]:
        color_by = trace.get("color_by")
        if color_by is not None:
            facet_and_color_by.add(color_by)

    if len(facet_and_color_by) > 0:
        agg_expr = "+ ".join(
            [f"approx_count_distinct({col})" for col in facet_and_color_by]
        )

        res = conn.table(table_name).aggregate(agg_expr).fetchone()

        if res is not None and res[0] > 50:
            # todo(rteqs): send error to frontend
            raise ValueError("Too many facets")

    relations: list[duckdb.DuckDBPyRelation] = []

    for i, trace in enumerate(config["traces"]):
        x = trace["x"]
        y = trace["y"]
        color_by = trace.get("color_by")
        error_bar = trace.get("error_bar")

        distinct_cols = (
            f"{x}, {y}"
            + (f", {facet}" if facet is not None else "")
            + (f", {color_by}" if color_by is not None else "")
        )

        cols = (
            distinct_cols
            + (f", {error_bar}" if error_bar not in {None, "sem", "stddev"} else "")
            + (f", {marker_size_axis}" if marker_size_axis is not None else "")
            + (
                f", array_value({custom_data_str})"
                if custom_data_str is not None
                else ""
            )
        )

        # note: can't do parameterized arguments for column and table names
        trace_data = conn.sql(
            f"""
            select
                distinct on ({distinct_cols})
                {cols}
            from
                {table_name}
            """
        ).set_alias(f"trace_{i}")

        # todo(rteqs): handle categorical axis
        min_max = trace_data.aggregate(
            f"min({x}) as min_x, min({y}) as min_y, max({x}) as max_x, max({y}) as max_y"
        ).set_alias("min_max")

        cell_size = 2
        max_occupancy = 3

        trace_data = (
            trace_data.join(min_max, "1 = 1")
            .project(
                f"""
                *,
                row_number() over (
                    partition by
                        floor((({x} - min_x) / (max_x - min_x) * {width_px}) / {cell_size}) +
                        floor((({y} - min_y) / (max_y - min_y) * {height_px}) / {cell_size}) * ceil({width_px} / {cell_size})
                ) as row_num
            """
            )
            .filter(f"row_num <= {max_occupancy}")
        )

        # todo(rteqs): slow to join with very large number of points, but if you have that many groups on the x-axis, you probably aren't using error bars
        if error_bar == "sem":
            sem = (
                conn.table(table_name)
                .aggregate(f"{x}, stddev({y} / sqrt(count({y}))) as sem")
                .set_alias(f"sem_{i}")
            )
            trace_data = trace_data.join(sem, f"trace_{i}.{x} = sem_{i}.{x}")

        elif error_bar == "stddev":
            stddev = (
                conn.table(table_name)
                .aggregate(f"{x}, stddev({y}) as stddev")
                .set_alias(f"stddev_{i}")
            )
            trace_data = trace_data.join(stddev, f"trace_{i}.{x} = stddev_{i}.{x}")

        relations.append(trace_data)

    return relations


async def downsample_ldata(
    conn: DuckDBPyConnection, ldata_node_id: str, config: PlotConfig
) -> list[list[Any]]:
    is_latest, last_modified_time = await check_generation(
        conn, f"ldata_{ldata_node_id}", "ldata", None
    )

    if not is_latest:
        url = await get_presigned_url(f"latch://{ldata_node_id}.node")
        conn.read_csv(url).to_table(f"ldata_{ldata_node_id}")
        conn.execute(
            """
            insert into
                plots_faas_catalog
            values
                ($name, null, $last_modified_time)
            on conflict
                (name)
            do update
            set
                last_modified_time = $last_modified_time
            """,
            parameters={
                "name": f"ldata_{ldata_node_id}",
                "last_modified_time": last_modified_time,
            },
        )

    return [rel.fetchall() for rel in downsample(conn, ldata_node_id, config)]


async def downsample_fig(
    conn: DuckDBPyConnection, key: str, fig: BaseFigure, cur_gen: int
) -> BaseFigure:
    is_latest, _ = await check_generation(conn, key, cur_gen=cur_gen)
    if not is_latest:
        # todo(rteqs): figure out how to convert plotly fig to duckdb table
        # df
        # conn.register(key, df)
        conn.execute(
            """
            insert into
                plots_faas_catalog
            values
                ($name, $generation, null)
            on conflict
                (name)
            do update
            set
                generation = $generation
            """,
            parameters={"name": key, "generation": cur_gen},
        )

    # convert fig to downsample config
    config: PlotConfig = {"traces": [{"x": "x", "y": "y"}]}
    relations = downsample(conn, key, config)


async def downsample_df(
    conn: DuckDBPyConnection, key: str, df: DataFrame, config: PlotConfig, cur_gen: int
) -> list[list[Any]]:
    is_latest, _ = await check_generation(conn, key, cur_gen=cur_gen)
    if not is_latest:
        conn.register(key, df.assign(index=df.index))
        conn.execute(
            """
            insert into
                plots_faas_catalog
            values
                ($name, $generation, null)
            on conflict
                (name)
            do update
            set
                generation = $generation
            """,
            parameters={"name": key, "generation": cur_gen},
        )

    return [rel.fetchall() for rel in downsample(conn, key, config)]
