# todo(rteqs): refactor this & move to new infra

import asyncio
import math
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import NotRequired, TypedDict

import duckdb
import numpy as np
from latch_asgi.auth import get_signer_sub
from latch_asgi.context.websocket import Context, HandlerResult
from latch_asgi.framework.websocket import WebsocketBadMessage
from latch_postgres.postgres import LatchAsyncConnection, sqlq
from psycopg import sql

from ..app.db import with_conn_retry

min_points = 100_000


@dataclass(frozen=True)
class AuthMessage:
    token: str
    ws_id: str


@dataclass(frozen=True)
class DbRes:
    access_ok: bool


@dataclass(frozen=True)
class PlotTrace:
    x: str
    y: str
    custom_data: list[str] | None
    color_by: str | None
    error_bar: str | None


@dataclass(frozen=True)
class SubsampleMessage:
    ldata_node_url: str
    ldata_node_id: str
    traces: list[PlotTrace]
    facet: str | None
    xrange: list[float | int] | None
    yrange: list[float | int] | None
    height_px: float
    width_px: float
    seq_num: int


class TraceData(TypedDict):
    x: list[list[float] | list[int] | list[object]]
    y: list[list[float] | list[int] | list[object]]
    custom_data: list[list[list[str]]]
    x_drank: NotRequired[list[list[int]]]
    y_drank: NotRequired[list[list[int]]]
    color_by: NotRequired[list[list[float] | list[int] | list[object]]]
    facet: NotRequired[list[list[float] | list[int] | list[object]]]
    error_bar: NotRequired[list[list[int | float]]]


@dataclass(frozen=True)
class PlotData:
    x: list[float] | list[int] | list[object]
    y: list[float] | list[int] | list[object]
    errorBar: list[int] | list[float] | None
    customData: list[list[str]]
    name: str
    xaxis: str
    yaxis: str


@dataclass(frozen=True)
class SubsampleResult:
    data: list[PlotData]
    subsampled: bool
    seqNum: int


class MaxGroupsExceededError(RuntimeError):
    """Raised when the number of groups exceeds the limit"""


class SubsampleError(RuntimeError):
    pass


def downsample(
    x: list[int] | list[float] | list[object],
    y: list[int] | list[float] | list[object],
    x_drank: list[int] | None,
    y_drank: list[int] | None,
    custom_data: list[object],
    xrange: tuple[float, float] | None,
    yrange: tuple[float, float] | None,
    width_px: float,
    height_px: float,
    error_bar: list[int] | list[float] | None,
) -> tuple[np.ndarray, np.ndarray, list[object], np.ndarray | None]:
    if len(x) != len(y):
        err_msg = (
            f"length of x and y must be equal, got len(x)={len(x)} and len(y)={len(y)}"
        )
        print(err_msg)
        raise SubsampleError(err_msg)

    x_col_np = np.array(x)
    y_col_np = np.array(y)
    x_drank_np = np.array(x_drank) if x_drank is not None else None
    y_drank_np = np.array(y_drank) if y_drank is not None else None
    custom_data = np.array(custom_data, dtype=object)
    error_bar_np = np.array(error_bar) if error_bar is not None else None

    x_np = x_col_np if x_drank_np is None else x_drank_np
    y_np = y_col_np if y_drank_np is None else y_drank_np

    mask = np.zeros(len(x_np), dtype=bool)

    # todo(rteqs): remove when we support string ranges. the argument
    if xrange is None:
        xrange = (np.min(x_np), np.max(x_np))

    if yrange is None:
        yrange = (np.min(y_np), np.max(y_np))

    xmin, xmax = sorted(xrange)
    ymin, ymax = sorted(yrange)

    x_np_scaled = (x_np - xmin) / (xmax - xmin) * width_px
    y_np_scaled = (y_np - ymin) / (ymax - ymin) * height_px
    points = np.column_stack((x_np_scaled, y_np_scaled))

    cell_size = 2
    max_occupancy = 3
    num_cols = math.ceil(width_px / cell_size)

    hash_ids = (
        np.floor_divide(x_np_scaled, cell_size)
        + np.floor_divide(y_np_scaled, cell_size) * num_cols
    ).astype(int)

    sort_idx = np.argsort(hash_ids)
    _, counts = np.unique(hash_ids[sort_idx], return_counts=True)

    k = -max_occupancy
    start = 0

    for count in counts:
        end = start + count

        bucket_indices = sort_idx[start:end]
        if len(bucket_indices) <= max_occupancy:
            mask[bucket_indices] = True

        else:
            points_in_bucket = points[bucket_indices]
            distances = np.sum((points_in_bucket - points_in_bucket[0]) ** 2, axis=1)
            max_dist_indices = np.argpartition(distances, k)[k:]
            mask_indices = bucket_indices[max_dist_indices]
            mask[mask_indices] = True
        start += count

    return (
        x_col_np[mask],
        y_col_np[mask],
        custom_data[mask],
        error_bar_np[mask] if error_bar_np is not None else None,
    )


def query_data(msg: SubsampleMessage) -> list[TraceData]:
    facet = msg.facet
    ldata_node_url = sql.Literal(msg.ldata_node_url)
    ldata_node_id = msg.ldata_node_id
    traces = msg.traces
    table_name = f"table_{ldata_node_id}"
    table_name_identifier = sqlq(table_name)

    xrange = msg.xrange
    yrange = msg.yrange

    data_dir = Path(Path.cwd() / ".data")
    if not Path.exists(data_dir):
        Path.mkdir(data_dir)
    parquet_file = Path(data_dir / f"{table_name}.parquet")
    parquet_file_literal = sql.Literal(str(parquet_file))

    trace_data: list[TraceData] = []

    with duckdb.connect() as conn:
        # todo(rteqs): versioning & cache / evicting
        # race condition ??
        conn.execute("set enable_progress_bar = false")
        if not Path.exists(parquet_file):
            create_parquet_q = (
                sqlq(
                    "copy (select *, row_number() over () as rank from read_csv_auto({})) to {} (format parquet)"
                )
                .format(ldata_node_url, parquet_file_literal)
                .as_string(None)
            )
            conn.execute(create_parquet_q)

        # todo(rteqs): not ideal since we read from the local file to create the table on every query.
        # the latency here is mitigated since parquet reads are roughly 2 orders fo magnitude faster than csv reads
        # ideally we write to a single .db file and read from that but this requires some thinking on concurrency
        # since duckdb only supports a sigle read-write connection or multiple read only connections.
        create_table_q = (
            sqlq("create table {} as (select * from {})")
            .format(table_name_identifier, parquet_file_literal)
            .as_string(None)
        )
        conn.execute(create_table_q)

        table_info_q = (
            sqlq(
                "select column_name, data_type from information_schema.columns where table_name = {}"
            )
            .format(sql.Literal(table_name))
            .as_string(None)
        )
        table_info_raw = conn.execute(table_info_q).fetchall()
        table_info = defaultdict(
            lambda: "BIGINT", {col: data_type for col, data_type in table_info_raw}
        )

        facet_and_color_by = [
            trace.color_by for trace in traces if trace.color_by is not None
        ] + ([facet] if facet is not None else [])

        if len(facet_and_color_by) > 0:
            select_parts = sqlq("+ ").join(
                [
                    sqlq('approx_count_distinct("{}")').format(sqlq(col))
                    for col in facet_and_color_by
                ]
            )
            check_max_num_groups_q = (
                sqlq(
                    """
                    select
                        {}
                    from
                        {}
                """
                )
                .format(select_parts, table_name_identifier)
                .as_string(None)
            )

            res = conn.execute(check_max_num_groups_q).fetchall()

            if res[0][0] > 50:
                raise MaxGroupsExceededError(
                    "Number of facet / color by groups exceeded the limit of 50"
                )

        for trace in traces:
            if (trace.x is None) or (trace.y is None):
                continue

            x_col = "rank" if trace.x == "__rank__" else trace.x
            x_in = sqlq('"{}" as x').format(sqlq(x_col))
            y_in = sqlq('"{}" as y').format(sqlq(trace.y))
            x_drank_in = (
                sqlq("dense_rank() over (order by x) as x_drank")
                if table_info[x_col] not in ["BIGINT", "DOUBLE"]
                else None
            )
            y_drank_in = (
                sqlq("dense_rank() over (order by y) as y_drank")
                if table_info[trace.y] not in ["BIGINT", "DOUBLE"]
                else None
            )

            custom_data_agg = sqlq("array[] as custom_data")
            if trace.custom_data is not None:
                custom_data_agg_parts = [
                    sqlq('list_distinct(array_agg(cast("{}" as varchar)))').format(
                        sqlq(col)
                    )
                    for col in trace.custom_data
                ]
                custom_data_agg = sqlq("array[{}] as custom_data").format(
                    sqlq(", ").join(custom_data_agg_parts)
                )

            color_by_in = (
                sqlq('"{}" as color_by').format(sqlq(trace.color_by))
                if trace.color_by is not None
                else None
            )

            error_bar_in = (
                sqlq('"{}" as error_bar').format(sqlq(trace.error_bar))
                if trace.error_bar is not None
                and table_info[trace.error_bar] in ["DOUBLE", "BIGINT"]
                else None
            )

            facet_in = (
                sqlq("{} as facet").format(sqlq(facet)) if facet is not None else None
            )

            inner_selection = sqlq(", ").join(
                [
                    col
                    for col in [
                        x_in,
                        y_in,
                        custom_data_agg,
                        color_by_in,
                        error_bar_in,
                        facet_in,
                        x_drank_in,
                        y_drank_in,
                    ]
                    if col is not None
                ]
            )

            inner_group_by_parts = [sqlq("x"), sqlq("y")]
            if trace.color_by is not None:
                inner_group_by_parts.append(sqlq("color_by"))
            if facet is not None:
                inner_group_by_parts.append(sqlq("facet"))
            if trace.error_bar is not None:
                inner_group_by_parts.append(sqlq("error_bar"))

            inner_group_by = sqlq("group by {}").format(
                sqlq(", ").join(inner_group_by_parts)
            )

            inner_where_parts = []
            # todo(rteqs): support string ranges
            if xrange is not None and table_info[x_col] in ["BIGINT", "DOUBLE"]:
                inner_where_parts.append(
                    sqlq("x >= {} and x <= {}").format(xrange[0], xrange[1])
                )
            if yrange is not None and table_info[trace.y] in ["BIGINT", "DOUBLE"]:
                inner_where_parts.append(
                    sqlq("y >= {} and y <= {}").format(yrange[0], yrange[1])
                )

            inner_where = (
                sqlq("and {}").format(sqlq(" and ").join(inner_where_parts))
                if len(inner_where_parts) > 0
                else sqlq("")
            )

            order_outer_selections_parts = [
                sqlq(col_name)
                for col_name, col in {
                    "x_drank": x_drank_in,
                    "y_drank": y_drank_in,
                }.items()
                if col is not None
            ]
            order_outer_selections = (
                sqlq("order by {}").format(
                    sqlq(", ").join(order_outer_selections_parts)
                )
                if len(order_outer_selections_parts) > 0
                else sqlq("")
            )

            outer_selection = sqlq(
                """
                array_agg(x {}) as x,
                array_agg(y {}) as y,
                array_agg(custom_data {}) as custom_data,
                {}
                {}
                {}
                {}
                {}
            """
            ).format(
                order_outer_selections,
                order_outer_selections,
                order_outer_selections,
                (
                    sqlq("array_agg(error_bar {}) as error_bar,").format(
                        order_outer_selections
                    )
                    if error_bar_in is not None
                    else sqlq("")
                ),
                (
                    sqlq("array_agg(x_drank {}) as x_drank,").format(
                        order_outer_selections
                    )
                    if x_drank_in is not None
                    else sqlq("")
                ),
                (
                    sqlq("array_agg(y_drank {}) as y_drank,").format(
                        order_outer_selections
                    )
                    if y_drank_in is not None
                    else sqlq("")
                ),
                sqlq("color_by,") if trace.color_by is not None else sqlq(""),
                sqlq("facet") if facet is not None else sqlq(""),
            )

            outer_group_by_parts = [
                sqlq(col_name)
                for col_name, col in {
                    "color_by": trace.color_by,
                    "facet": facet,
                }.items()
                if col is not None
            ]
            outer_group_by = (
                sqlq('group by "{}"').format(sqlq(", ").join(outer_group_by_parts))
                if len(outer_group_by_parts) > 0
                else sqlq("")
            )

            order_by_parts = [
                sqlq(col_name)
                for col_name, col in {
                    "color_by": trace.color_by,
                    "facet": facet,
                }.items()
                if col is not None
            ]
            order_by = (
                sqlq('order by "{}"').format(sqlq(", ").join(order_by_parts))
                if len(order_by_parts) > 0
                else sqlq("")
            )

            query = (
                sqlq(
                    """
                    select
                        {}
                    from
                        (
                            select
                                {}
                            from
                                {}
                            where
                                x is not null and
                                y is not null
                                {}
                            {}
                        )
                    {}
                    {}
                """
                )
                .format(
                    outer_selection,
                    inner_selection,
                    table_name_identifier,
                    inner_where,
                    inner_group_by,
                    outer_group_by,
                    order_by,
                )
                .as_string(None)
            )

            res = conn.execute(query).fetchnumpy()
            trace_data.append(res)

    return trace_data


async def get_latest_message(ctx: Context) -> SubsampleMessage:
    msg = await ctx.receive_message(SubsampleMessage)

    # if we've recieved more then one message, we get the latest one and drop the rest
    ten_ms = 1 / 100
    while True:
        try:
            msg = await asyncio.wait_for(
                ctx.receive_message(SubsampleMessage), timeout=ten_ms
            )
        except TimeoutError:
            break

    return msg


def run_task(msg: SubsampleMessage) -> SubsampleResult:
    facet_groups: list[str] = []
    plot_data: list[PlotData] = []
    subsampled = False

    data = query_data(msg)

    for trace in data:
        for i in range(len(trace["x"])):
            x = trace["x"][i]
            y = trace["y"][i]
            x_drank = trace["x_drank"][i] if "x_drank" in trace else None
            y_drank = trace["y_drank"][i] if "y_drank" in trace else None
            custom_data = trace["custom_data"][i]
            color_by = trace["color_by"][i] if "color_by" in trace else None
            facet = trace["facet"][i] if "facet" in trace else None
            error_bar = trace["error_bar"][i] if "error_bar" in trace else None
            if facet is not None:
                if not isinstance(facet, str):
                    facet = facet.item()
                facet_groups.append(facet)

            name = f"Series {i}"
            if facet is not None:
                name = str(facet)
            elif color_by is not None:
                name = str(color_by)

            axis_idx = "" if i == 0 else i + 1

            if len(x) <= min_points:
                plot_data.append(
                    {
                        "x": x,
                        "y": y,
                        "customData": custom_data,
                        "name": name,
                        "errorBar": error_bar,
                        "xaxis": (f"x{axis_idx}" if facet is not None else None),
                        "yaxis": (f"y{axis_idx}" if facet is not None else None),
                    }
                )
                continue

            subsampled = True
            subsampled_data = downsample(
                x,
                y,
                x_drank,
                y_drank,
                custom_data,
                msg.xrange,
                msg.yrange,
                msg.width_px,
                msg.height_px,
                error_bar,
            )

            plot_data.append(
                {
                    "x": subsampled_data[0].tolist(),
                    "y": subsampled_data[1].tolist(),
                    "customData": subsampled_data[2].tolist(),
                    "errorBar": (
                        subsampled_data[3].tolist()
                        if subsampled_data[3] is not None
                        else None
                    ),
                    "name": name,
                    "xaxis": (f"x{axis_idx}" if facet is not None else None),
                    "yaxis": (f"y{axis_idx}" if facet is not None else None),
                }
            )

    return {
        "data": {"plotData": plot_data, "facetGroups": facet_groups},
        "subsampled": subsampled,
        "seqNum": msg.seq_num,
    }


executor = ProcessPoolExecutor()


async def subsample(ctx: Context) -> HandlerResult:
    await ctx.accept_connection()

    auth_msg = await ctx.receive_message(AuthMessage)
    auth = get_signer_sub(auth_msg.token)

    @with_conn_retry
    async def db_work(conn: LatchAsyncConnection) -> DbRes:
        async with conn.transaction():
            return await conn.query1(
                DbRes,
                sqlq(
                    """
                    select
                        app_private.is_account_viewer_non_authenticated(
                            %(ws_id)s,
                            sub_acc_id,
                            'member'
                        )
                            as access_ok
                    from
                        app_private.backend_resolve_auth(
                            %(sub)s,
                            %(exec_tok)s,
                            %(sdk_tok)s
                        )
                        _(sub_acc_id)
                    """
                ),
                sub=auth.oauth_sub,
                exec_tok=auth.execution_token,
                sdk_tok=auth.sdk_token,
                ws_id=auth_msg.ws_id,
            )

    db_res = await db_work()
    if not db_res.access_ok:
        raise WebsocketBadMessage(f"Signer cannot access workspace {auth_msg.ws_id}")

    loop = asyncio.get_event_loop()
    while True:
        msg = await get_latest_message(ctx)

        try:
            res = await loop.run_in_executor(executor, run_task, msg)
            await ctx.send_message(res)

        except MaxGroupsExceededError as e:
            await ctx.send_message({"error": str(e), "seqNum": msg.seq_num})

        except SubsampleError as e:
            # todo(rteqs): attach to span
            print("error:", e)
            await ctx.send_message({"error": "Subsample Error", "seqNum": msg.seq_num})
            break
