from time import time
from typing import Literal, NotRequired, TypedDict

import duckdb
from lplots import _inject
from pandas import DataFrame
from plotly.basedatatypes import BaseFigure

conn = duckdb.connect(
    database=":memory:plots-faas",
    read_only=False,
    config={"enable_progress_bar": False},
)

conn.execute(
    """
    create table plots_faas_catalog (
        name string not null unique,
        generation integer,
        last_modified_time timestamp
    )
    """
)


class Trace(TypedDict):
    x: str
    y: str
    color_by: NotRequired[str]
    error_bar: NotRequired[str]


class DownsamplePlotConfig(TypedDict):
    traces: list[Trace]
    custom_data: NotRequired[list[str]]

    facet: NotRequired[str]
    marker_size_axis: NotRequired[str]

    xrange: NotRequired[tuple[float | int, float | int]]
    yrange: NotRequired[tuple[float | int, float | int]]
    height_px: NotRequired[float]
    width_px: NotRequired[float]


def has_latest_generation(
    table_name: str, source_type: Literal["kernel", "ldata"] = "kernel"
) -> bool:
    cur_gen = (
        _inject.kernel.k_globals.generation_counter.get(table_name)
        if source_type == "kernel"
        else None
    )

    # todo(rteqs): ldata_last_modified_time gql query to check latest ldata_node_id event
    last_modified_time = time()

    duckdb_gen = conn.sql(
        """
        select
            generation = $cur_gen or
            last_modified_time > $last_modified_time
        from
            plots_faas_catalog
        where
            name = $table_name
        """,
        params={
            "table_name": table_name,
            "cur_gen": cur_gen,
            "last_modified_time": last_modified_time,
        },
    ).fetchone()

    if duckdb_gen is None or not isinstance(duckdb_gen[0], bool):
        raise ValueError(f"Unexpected generation value: {duckdb_gen}")

    return duckdb_gen[0]


def downsample(
    table_name: str, config: DownsamplePlotConfig
) -> list[duckdb.DuckDBPyRelation]:
    custom_data = config.get("custom_data")
    custom_data_str = ", ".join(custom_data) if custom_data is not None else None

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
        select_parts = "+ ".join(
            [f"approx_count_distinct({col})" for col in facet_and_color_by]
        )

        res = conn.sql(
            """
            select
                $select_parts
            from
                $table_name
            """,
            params={"table_name": table_name, "select_parts": select_parts},
        ).fetchone()

        if res is not None and res[0] > 50:
            # todo(rteqs): wip
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

        trace_data = conn.sql(
            """
            select
                distinct on ($distinct_cols)
                $cols
            from
                $table_name
            """,
            params={
                "table_name": table_name,
                "distinct_cols": distinct_cols,
                "cols": cols,
            },
        ).set_alias(f"trace_{i}")

        # todo(rteqs): handle categorical axis
        min_max = trace_data.project(
            f"min({x}) as min_x, min({y}) as min_y, max({x}) as max_x, max({y}) as max_y"
        )

        cell_size = 2
        max_occupancy = 3

        trace_data = (
            trace_data.join(min_max, condition="1 = 1")  # cross join
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
            sem = conn.sql(
                """
                select
                    $x,
                    stddev($y) / sqrt(count($y)) as sem
                from
                    $table_name
                group by
                    $x
                """,
                params={"table_name": table_name, "x": x, "y": y},
            ).set_alias(f"sem_{i}")
            trace_data = trace_data.join(sem, f"trace_{i}.{x} = sem_{i}.{x}")

        elif error_bar == "stddev":
            stddev = conn.sql(
                """
                select
                    $x,
                    stddev($y) as stddev
                from
                    $table_name
                group by
                    $x
                """,
                params={"table_name": table_name, "x": x, "y": y},
            ).set_alias("stddev")
            trace_data = trace_data.join(stddev, f"trace_{i}.{x} = stddev_{i}.{x}")

        relations.append(trace_data)

    return relations


def downsample_ldata(
    ldata_node_id: str, config: DownsamplePlotConfig
) -> list[DataFrame]:
    if not has_latest_generation(f"ldata_{ldata_node_id}", "ldata"):
        pass

    # todo(rteqs): process json on our own to avoid extra cost of going through pandas
    return [rel.df() for rel in downsample(ldata_node_id, config)]


# def downsample(key: str):


def downsample_fig(key: str, fig: BaseFigure) -> BaseFigure:
    if not has_latest_generation(key):
        # todo(rteqs): import
        pass

    # convert fig to downsample config
    config: DownsamplePlotConfig = {"traces": [{"x": "x", "y": "y"}]}
    relations = downsample(key, config)


def downsample_df(
    key: str, df: DataFrame, config: DownsamplePlotConfig
) -> list[DataFrame]:
    if not has_latest_generation(key):
        pass

    # todo(rteqs): process json on our own to avoid extra cost of going through pandas
    return [rel.df() for rel in downsample(key, config)]
