import datetime
import sys
import traceback
from pathlib import Path
from typing import Any, Literal, TypedDict

from duckdb import ColumnExpression as Col
from duckdb import ConstantExpression as Const
from duckdb import DuckDBPyConnection, DuckDBPyRelation
from duckdb import connect as duckdb_connect
from pandas import DataFrame
from plotly.basedatatypes import BaseFigure

# todo(rteqs): get rid of this
sys.path.append(str(Path(__file__).parent.absolute()))
from utils import PlotConfig, auth_token_sdk, get_presigned_url, gql_query


def initialize_duckdb() -> DuckDBPyConnection:
    conn = duckdb_connect(database=":memory:plots-faas", read_only=False)

    conn.execute(
        """
        create table plots_faas_catalog (
            name string not null unique,
            in_memory_data_generation integer,
            external_data_last_modified_time timestamp
            check (in_memory_data_generation is not null or external_data_last_modified_time is not null)
        )
        """
    )
    return conn


# todo(rteqs): replace with psycopg quote. this would require base image update to install psycopg
def quote_identifier(x: str) -> str:
    assert '"' not in x
    return f'"{x}"'


class DownsampleResult(TypedDict):
    columns: list[Any]
    data: list[Any]


async def get_staleness_info(
    conn: DuckDBPyConnection,
    table_name: str,
    source_type: Literal["kernel", "ldata"] = "kernel",
    cur_gen: int | None = None,
) -> tuple[bool, datetime.datetime | None]:
    external_data_last_modified_time = None
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
                external_data_last_modified_time = datetime.datetime.fromisoformat(
                    data[0]["time"]
                )
        except Exception:
            traceback.print_exc()
            raise

    external_data_last_modified_time = datetime.datetime.now(tz=datetime.UTC)

    # fixme(rteqs): whut? transaction aborted error. might need some kind of lock
    retries = 0
    max_retries = 10
    duckdb_gen = None

    while retries < max_retries:
        try:
            duckdb_gen = (
                conn.table("plots_faas_catalog")
                .filter(Col("name") == Const(table_name))
                .project(
                    (Col("in_memory_data_generation") == Const(cur_gen))
                    or (
                        Col("external_data_last_modified_time")
                        >= Const(external_data_last_modified_time)
                    )
                )
                .fetchone()
            )
            break
        except Exception:
            conn.rollback()
            retries += 1
            if retries == max_retries:
                raise

    if duckdb_gen is None:
        return False, None

    if not isinstance(duckdb_gen[0], bool):
        raise TypeError(f"Unexpected generation value: {duckdb_gen}")

    return duckdb_gen[0], external_data_last_modified_time


def downsample(
    conn: DuckDBPyConnection, table_name: str, config: PlotConfig
) -> list[DuckDBPyRelation]:
    custom_data = [quote_identifier(col) for col in config.get("custom_data", [])]
    custom_data_str = ", ".join(custom_data) if len(custom_data) > 0 else None

    facet = config.get("facet")
    facet = quote_identifier(facet) if facet is not None else None

    if facet is not None:
        res = (
            conn.table(table_name)
            .aggregate(f"approx_count_distinct({facet})")
            .fetchone()
        )

        if res is not None and res[0] > 50:
            # todo(rteqs): send error to frontend
            raise ValueError("Too many facets")

    width_px = config.get("width_px", 2000)
    height_px = config.get("height_px", 600)

    [min_x, max_x] = config.get("xrange", [None, None])
    [min_y, max_y] = config.get("yrange", [None, None])

    relations: list[DuckDBPyRelation] = []

    categorical_columns = set(
        conn.query(
            """
            select
                array_agg(column_name)
            from
                information_schema.columns
            where
                table_name=$table_name and
                data_type not in ['DOUBLE', 'BIGINT']
            """,
            params={"table_name": table_name},
        ).fetchall()[0][0]
    )

    for i, trace in enumerate(config["traces"]):
        if trace["type"] != "scattergl" and trace["type"] != "scatter":
            continue

        x_raw = trace["x"]
        y_raw = trace["y"]

        x_original = x = quote_identifier(trace["x"])
        y_original = y = quote_identifier(trace["y"])

        x_is_categorical = x_raw in categorical_columns
        y_is_categorical = y_raw in categorical_columns

        color_by = trace.get("color_by")
        color_by = quote_identifier(color_by) if color_by is not None else None

        error_bar = trace.get("error_bar")
        error_bar = quote_identifier(error_bar) if error_bar is not None else None

        marker_size = trace.get("marker_size")
        marker_size = quote_identifier(marker_size) if marker_size is not None else None

        distinct_cols = (
            f"{x}, {y}"
            + (f", {facet}" if facet is not None else "")
            + (f", {color_by}" if color_by is not None else "")
        )

        cols = (
            distinct_cols
            + (f", {min_x} as min_x" if min_x is not None else "")
            + (f", {max_x} as max_x" if max_x is not None else "")
            + (f", {min_y} as min_y" if min_y is not None else "")
            + (f", {max_y} as max_y" if max_y is not None else "")
            + (f", {error_bar}" if error_bar not in {None, "sem", "stddev"} else "")
            + (f", {marker_size}" if marker_size is not None else "")
            + (
                f", array_value({custom_data_str}) as custom_data"
                if custom_data_str is not None
                else ""
            )
            + (
                f", dense_rank() over (order by {x}) as x_dense_rank"
                if x_is_categorical
                else ""
            )
            + (
                f", dense_rank() over (order by {y}) as y_dense_rank"
                if y_is_categorical
                else ""
            )
        )

        x = "x_dense_rank" if x_is_categorical else x
        y = "y_dense_rank" if y_is_categorical else y

        trace_data = conn.sql(
            f"""
            select
                distinct on ({distinct_cols})
                {cols}
            from
                {table_name}
            """
        ).set_alias(f"trace_{i}")

        agg_expr = ""
        for col, val in [
            (f"min({x}) as min_x", min_x),
            (f"max({x}) as max_x", max_x),
            (f"min({y}) as min_y", min_y),
            (f"max({y}) as max_y", max_y),
            # extrema
            (f"min({x}) as global_min_x", None),
            (f"max({x}) as global_max_x", None),
            (f"min({y}) as global_min_y", None),
            (f"max({y}) as global_max_y", None),
        ]:
            if val is None:
                agg_expr += f", {col}" if len(agg_expr) > 0 else col

        excluded_columns = [
            "row_num",
            "min_x",
            "max_x",
            "min_y",
            "max_y",
            "global_min_x",
            "global_max_x",
            "global_min_y",
            "global_max_y",
        ]

        if x_is_categorical:
            excluded_columns.append("x_dense_rank")
        if y_is_categorical:
            excluded_columns.append("y_dense_rank")

        exclude_str = ", ".join(excluded_columns)

        cell_size = 3
        max_occupancy = 2

        min_max = trace_data.aggregate(agg_expr).set_alias("min_max")
        trace_data = (
            trace_data.join(min_max, "1 = 1")
            .filter(
                f"""
                    (
                        min_x <= {x} and {x} <= max_x and
                        min_y <= {y} and {y} <= max_y
                    )
                    or {x} = global_min_x
                    or {x} = global_max_x
                    or {y} = global_min_y
                    or {y} = global_max_y
                """
            )
            .project(
                f"""
                *,
                row_number() over (
                    partition by
                        {f"{facet}, " if facet is not None else ""}
                        floor((({x} - min_x) / (max_x - min_x) * {width_px}) / {cell_size}) +
                        floor((({y} - min_y) / (max_y - min_y) * {height_px}) / {cell_size}) * ceil({width_px} / {cell_size})
                        {f", {color_by}" if color_by is not None else ""}
                ) as row_num
            """
            )
            .filter(f"row_num <= {max_occupancy}")
            .project(f"* exclude({exclude_str})")
            .order(
                ",".join(
                    col
                    for col in [
                        facet,
                        color_by,
                        x_original if x_is_categorical else None,
                        y_original if y_is_categorical else None,
                    ]
                    if col is not None
                )
            )
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
) -> list[DownsampleResult]:
    # todo(rteqs): wip
    is_latest, external_data_last_modified_time = await get_staleness_info(
        conn, f"ldata_{ldata_node_id}", "ldata", None
    )

    if not is_latest:
        url = await get_presigned_url(f"latch://{ldata_node_id}.node")
        conn.read_csv(url).to_table(f"ldata_{ldata_node_id}")
        conn.execute(
            """
            insert into
                plots_faas_catalog (
                    name,
                    external_data_last_modified_time
                )
            values
                ($name, $external_data_last_modified_time)
            on conflict
                (name)
            do update
            set
                external_data_last_modified_time = $external_data_last_modified_time
            """,
            parameters={
                "name": f"ldata_{ldata_node_id}",
                "external_data_last_modified_time": external_data_last_modified_time,
            },
        )

    return [
        {"columns": rel.columns, "data": rel.fetchall()}
        for rel in downsample(conn, f"ldata_{ldata_node_id}", config)
    ]


async def downsample_fig(
    conn: DuckDBPyConnection, key: str, fig: BaseFigure, cur_gen: int
) -> BaseFigure:
    table_name = f"in_memory_{key}"
    is_latest, _ = await get_staleness_info(conn, table_name, cur_gen=cur_gen)
    if not is_latest:
        # todo(rteqs): figure out how to convert plotly fig to duckdb table
        # df
        # conn.register(key, df)
        conn.execute(
            """
            insert into
                plots_faas_catalog (name, in_memory_data_generation)
            values
                ($name, $in_memory_data_generation)
            on conflict
                (name)
            do update
            set
                in_memory_data_generation = $in_memory_data_generation
            """,
            parameters={"name": table_name, "in_memory_data_generation": cur_gen},
        )

    # todo(rteqs): convert fig to downsample config
    config: PlotConfig = {"traces": [{"x": "x", "y": "y"}]}
    relations = downsample(conn, table_name, config)


async def downsample_df(
    conn: DuckDBPyConnection, key: str, df: DataFrame, config: PlotConfig, cur_gen: int
) -> list[DownsampleResult]:
    table_name = f"in_memory_{key}"
    is_latest, _ = await get_staleness_info(conn, table_name, cur_gen=cur_gen)
    if not is_latest:
        conn.register(table_name, df.assign(index=df.index))

        retries = 0
        max_retries = 10

        while retries < max_retries:
            try:
                conn.begin()
                conn.execute(
                    """
                    insert into
                        plots_faas_catalog (name, in_memory_data_generation)
                    values
                        ($name, $in_memory_data_generation)
                    on conflict
                        (name)
                    do update
                    set
                        in_memory_data_generation = $in_memory_data_generation
                    """,
                    parameters={
                        "name": table_name,
                        "in_memory_data_generation": cur_gen,
                    },
                )
                conn.commit()
                break
            except Exception:
                conn.rollback()
                retries += 1
                if retries == max_retries:
                    raise

    return [
        {"columns": rel.columns, "data": rel.fetchall()}
        for rel in downsample(conn, table_name, config)
    ]
