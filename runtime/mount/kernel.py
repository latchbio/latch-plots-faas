import ast
import asyncio
import math
import pprint
import re
import signal
import socket
import sys
import traceback
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from io import TextIOWrapper
from pathlib import Path
from traceback import format_exc
from types import FrameType
from typing import TYPE_CHECKING, Any, Literal, TypedDict, TypeVar

import numpy as np
import orjson
import pandas as pd
import plotly.io as pio
import plotly.io._json as pio_json
from duckdb import DuckDBPyConnection
from latch.registry.table import Table
from lplots import _inject
from lplots.reactive import Node, Signal, ctx
from lplots.themes import graphpad_inspired_theme
from lplots.utils.nothing import Nothing
from lplots.widgets._emit import WidgetState
from matplotlib.figure import Figure
from pandas import DataFrame, Index, MultiIndex, Series
from pandas.io.json._table_schema import build_table_schema
from plotly.basedatatypes import BaseFigure
from plotly_utils.precalc_box import precalc_box
from plotly_utils.precalc_violin import precalc_violin
from socketio_thread import SocketIoThread
from stdio_over_socket import SocketWriter, text_socket_writer

sys.path.append(str(Path(__file__).parent.absolute()))
from subsample import downsample_df, initialize_duckdb, quote_identifier
from utils import PlotConfig, get_presigned_url

if TYPE_CHECKING:
    from numpy.typing import NDArray

sys.path.pop()

assert isinstance(sys.stdout, TextIOWrapper)
sys.stdout.reconfigure(line_buffering=True)

assert isinstance(sys.stderr, TextIOWrapper)
sys.stderr.reconfigure(line_buffering=True)

K = TypeVar("K")
V = TypeVar("V")


class SortSetting(TypedDict):
    columnId: str
    direction: Literal["asc", "desc"]


class StringFilter(TypedDict):
    type: Literal["string"]
    value: str


class NumberFilter(TypedDict):
    type: Literal["number"]
    value: int | float


class BooleanFilter(TypedDict):
    type: Literal["boolean"]
    value: bool


class DatetimeFilter(TypedDict):
    type: Literal["datetime"]
    value: int


class ArrayFilter(TypedDict):
    type: Literal["array"]
    value: list["FilterValue"]


FilterValue = StringFilter | NumberFilter | BooleanFilter | DatetimeFilter | ArrayFilter

FilterOpCode = Literal[
    "empty",
    "=",
    ">=",
    "<=",
    "regex-contains",
    "regex-starts-with",
    "regex-ends-with",
    "one-of",
]


class FilterOperatorOptions(TypedDict):
    caseInsensitive: bool
    negate: bool


class FilterOperator(TypedDict):
    opcode: FilterOpCode
    options: FilterOperatorOptions


class Filter(TypedDict):
    columnId: str
    value: FilterValue
    operator: FilterOperator


ColSelections = tuple[str, list[str | int | float]]
DataframeSelections = list[tuple[ColSelections, ColSelections]]


@dataclass
class PaginationSettings:
    page_size: int = 25
    page_idx: int = 0
    sort_settings: SortSetting | None = None
    row_filters: list[Filter] | None = None
    selections: DataframeSelections | None = None


class TracedDict(dict[str, Signal[object]]):
    touched: set[str]
    removed: set[str]

    dataframes: Signal[set[str]]

    item_write_counter: defaultdict[str, int]
    duckdb: DuckDBPyConnection

    def __init__(self, duckdb: DuckDBPyConnection) -> None:
        self.touched = set()
        self.removed = set()

        self.dataframes = Signal(set())
        self.item_write_counter = defaultdict(int)
        self.duckdb = duckdb

    def __getitem__(self, __key: str) -> object:
        if __key == "__builtins__":
            return super().__getitem__("__builtins__")

        return self.getitem_signal(__key).sample()

    def getitem_signal(self, __key: str) -> Signal[object]:
        return super().__getitem__(__key)

    def get_signal(self, __key: str) -> Signal[object] | None:
        if __key not in self:
            return None

        return self.getitem_signal(__key)

    def __setitem__(self, __key: str, __value: object) -> None:
        self.touched.add(__key)
        self.item_write_counter[__key] += 1

        if __key == "__builtins__":
            return super().__setitem__(__key, __value)

        dfs = self.dataframes.sample()
        if hasattr(__value, "iloc") and __key not in dfs:
            dfs.add(__key)
            self.dataframes(dfs)

        if __key in self:
            sig = super().__getitem__(__key)

            old = sig.sample()
            if isinstance(__value, Signal) and isinstance(old, Signal):
                # allow simply setting `sig = Signal(val)`
                # without having to check `if "sig" in globals()`
                # to define a signal without losing its subscribers on re-runs
                old(__value.sample())
                return None

            sig(__value)
            sig._apply_updates()
            return None

        return super().__setitem__(__key, Signal(__value))

    def __delitem__(self, __key: str) -> None:
        self.touched.add(__key)
        self.removed.add(__key)
        if __key in self.item_write_counter:
            del self.item_write_counter[__key]
            table_name = f"in_memory_{quote_identifier(__key)}"
            self.duckdb.execute(
                f"""
                delete from
                    plots_faas_catalog
                where
                    name = {table_name}
                """
            )
            self.duckdb.execute(f"drop table {table_name}")

        dfs = self.dataframes.sample()
        if __key in dfs:
            dfs.remove(__key)
            self.dataframes(dfs)

        return super().__delitem__(__key)

    def clear(self) -> None:
        self.touched.clear()
        self.removed.clear()
        self.item_write_counter.clear()

    @property
    def available(self) -> set[str]:
        return self.touched - self.removed


class ExitException(Exception): ...


class CellInterruptException(Exception): ...


KeyType = Literal["key", "ldata_node_id", "registry_table_id", "url"]


def cell_exit(code: int = 0) -> None:
    raise ExitException


def cell_interrupt(code: int = 0) -> None:
    raise CellInterruptException


leading_digits_and_dash = re.compile(r"^\d+-")


def remove_leading_digits_and_dash(col_raw: str) -> str:
    return re.sub(leading_digits_and_dash, "", col_raw)


multi_index_col_name = re.compile(r"^level_\d+$")


def is_multi_index_col(col: str) -> bool:
    return re.match(multi_index_col_name, col) is not None


def filter_dataframe(
    df: DataFrame, col: str, op: FilterOperator, filter: FilterValue
) -> DataFrame:
    options: FilterOperatorOptions = op.get("options", {})
    case_insensitive = options.get("case_insensitive", False)
    negate = options.get("negate", False)

    opcode = op.get("opcode")
    filter_type = filter.get("type")
    filter_value = filter.get("value")

    col_vals: Series[Any] | Index[Any] | None = None
    if col == "index":
        col_vals = df.index
    elif is_multi_index_col(col) and isinstance(df.index, MultiIndex):
        level = int(col.split("_")[-1])
        col_vals = df.index.get_level_values(level)
    elif col in df.index.names:
        col_vals = df.index.get_level_values(col)
    else:
        if col not in df:
            return df

        col_vals = df[col]

    if opcode == "empty":
        return df[col_vals.notna() & (col_vals != "")]

    if filter_value is None:
        return df

    mask: Series[bool] | NDArray[np.bool] | None = None

    if opcode == "=":
        if np.issubdtype(col_vals.dtype.type, np.floating):
            mask = np.isclose(col_vals, np.array(filter_value))
        else:
            mask = col_vals == filter_value

    elif filter_type in {"number", "datetime"}:
        if opcode == ">=":
            mask = col_vals >= filter_value

        elif opcode == "<=":
            mask = col_vals <= filter_value

    elif filter_type == "string":
        assert isinstance(filter_value, str)

        # todo(rteqs): remove casting after we clean up frontend's mess
        if opcode == "regex-contains":
            mask = col_vals.astype(str).str.contains(
                filter_value, case=case_insensitive
            )

        elif opcode == "regex-starts-with":
            mask = col_vals.astype(str).str.startswith(filter_value)

        elif opcode == "regex-ends-with":
            mask = col_vals.astype(str).str.endswith(filter_value)

    elif filter_type == "array" and opcode == "one-of":
        assert isinstance(filter_value, list)

        if np.issubdtype(col_vals.dtype.type, np.floating):
            mask = np.zeros(len(df), dtype=bool)
            for val in filter_value:
                mask |= np.isclose(col_vals, np.array(val.get("value")))

        else:
            filter_value = [val.get("value") for val in filter_value]
            mask = col_vals.isin(filter_value)

    if negate and mask is not None:
        mask = ~mask

    res = df[mask]
    if TYPE_CHECKING:
        assert isinstance(res, DataFrame)

    return res


def filter_dataframe_by_selections(
    df: DataFrame, col: str, selections: list
) -> "DataFrame | Series[bool] | NDArray[np.bool]":
    if col not in df:
        return np.ones(len(df), dtype=bool)

    if np.issubdtype(df[col].dtype.type, np.floating):
        mask = np.zeros(len(df), dtype=bool)
        for sel in selections:
            mask |= np.isclose(df[col], sel)
        return mask

    return df[col].isin(selections)


def filter_and_sort(
    *, df: DataFrame, pagination_settings: PaginationSettings
) -> DataFrame:
    row_filters = pagination_settings.row_filters
    sort_settings = pagination_settings.sort_settings
    selections = pagination_settings.selections

    if row_filters is not None:
        for rf in row_filters:
            col_raw = rf.get("columnId")
            col = remove_leading_digits_and_dash(col_raw)
            op = rf.get("operator")
            value = rf.get("value")
            df = filter_dataframe(df, col, op, value)

    if (
        sort_settings is not None
        and "direction" in sort_settings
        and "columnId" in sort_settings
    ):
        direction = sort_settings.get("direction")
        col_raw = sort_settings.get("columnId")

        is_asc = direction == "asc"
        col = remove_leading_digits_and_dash(col_raw)

        if col == "index":
            df = df.sort_index(ascending=is_asc)
        elif is_multi_index_col(col) and isinstance(df.index, MultiIndex):
            df = df.sort_index(level=int(col.split("_")[-1]), ascending=is_asc)
        elif col in df.index.names:
            df = df.sort_index(level=col, ascending=is_asc)
        elif col in df:
            df = df.sort_values(by=col, ascending=is_asc)

    if selections is not None:
        if len(selections) == 0:
            df = pd.DataFrame()

        mask = np.zeros(len(df), dtype=bool)
        for trace in selections:
            sub_mask = np.ones(len(df), dtype=bool)
            for col, sel in trace:
                if col == "index":
                    sub_mask &= df.index.isin(sel)
                elif is_multi_index_col(col) and isinstance(df.index, MultiIndex):
                    level = int(col.split("_")[-1])
                    sub_mask &= df.index.get_level_values(level).isin(sel)
                elif col in df.index.names:
                    sub_mask &= df.index.get_level_values(col).isin(sel)
                else:
                    sub_mask &= filter_dataframe_by_selections(df, col, sel)
            mask |= sub_mask
        df = df[mask]

    return df


def paginate(*, df: DataFrame, pagination_settings: PaginationSettings) -> DataFrame:
    page_size = pagination_settings.page_size
    page_idx = pagination_settings.page_idx

    num_pages = math.ceil(len(df) / page_size)

    # if requested page index is out of bounds after filtering
    page_idx = max(0, page_idx)
    page_idx = min(page_idx, num_pages - 1)

    y = page_idx * page_size
    h = page_size

    # todo(rteqs): handle wide dataframes
    x = 0
    w = 100

    if hasattr(df, "compute"):
        # to support Dask and similar datasets which cannot index on rows
        data = df.iloc[:, x : x + w].head(w)
    else:
        data = df.iloc[y : y + h, x : x + w]

    return data


def pagination_settings_dict_factory() -> (
    defaultdict[str, defaultdict[str, PaginationSettings]]
):
    return defaultdict(lambda: defaultdict(PaginationSettings))


class DfJsonSplitFormat(TypedDict):
    columns: list[str]
    index: list[Any] | None
    data: list[list[Any]]


@dataclass
class CategorizedCellOutputs:
    all: list[str] = field(default_factory=list)
    dfs: list[str] = field(default_factory=list)
    figures: list[str] = field(default_factory=list)
    static_figures: list[str] = field(default_factory=list)


def serialize_plotly_figure(x: BaseFigure) -> object:
    res = x.to_dict()

    for trace in res["data"]:
        try:
            if trace["type"] == "box":
                precalc_box(trace)
            elif trace["type"] == "violin":
                precalc_violin(trace)
            elif trace["type"] == "scatter":
                trace["type"] = "scattergl"
        except Exception:
            traceback.print_exc()

    modules = {
        "sage_all": pio_json.get_module("sage.all", should_load=False),
        "np": pio_json.get_module("numpy", should_load=False),
        "pd": pio_json.get_module("pandas", should_load=False),
        "image": pio_json.get_module("PIL.Image", should_load=False),
    }

    # note(maximsmol): plotly itself does a bunch of escaping to avoid XSS
    # when embedding directly into HTML. we never do that so we don't care
    return pio_json.clean_to_json_compatible(res, modules=modules)


@dataclass(kw_only=True)
class Kernel:
    conn: SocketIoThread

    cell_seq = 0
    cell_rnodes: dict[str, Node] = field(default_factory=dict)
    k_globals: TracedDict = field(init=False)
    cell_status: dict[str, str] = field(default_factory=dict)

    active_cell: str | None = None

    widget_signals: dict[str, Signal[Any]] = field(default_factory=dict)
    nodes_with_widgets: dict[int, Node] = field(default_factory=dict)

    cell_output_selections: dict[str, str] = field(default_factory=dict)
    viewer_cell_selections: dict[str, tuple[str, KeyType]] = field(default_factory=dict)
    plot_data_selections: dict[str, str] = field(default_factory=dict)

    ldata_dataframes: dict[str, DataFrame] = field(default_factory=dict)
    registry_dataframes: dict[str, DataFrame] = field(default_factory=dict)
    url_dataframes: dict[str, DataFrame] = field(default_factory=dict)

    cell_pagination_settings: defaultdict[str, defaultdict[str, PaginationSettings]] = (
        field(default_factory=pagination_settings_dict_factory)
    )
    viewer_pagination_settings: defaultdict[
        str, defaultdict[str, PaginationSettings]
    ] = field(default_factory=pagination_settings_dict_factory)

    plot_configs: dict[str, PlotConfig | None] = field(default_factory=dict)
    duckdb: DuckDBPyConnection = field(default=initialize_duckdb())

    def __post_init__(self) -> None:
        self.k_globals = TracedDict(self.duckdb)
        self.k_globals["exit"] = cell_exit
        self.k_globals.clear()
        pio.templates["graphpad_inspired_theme"] = graphpad_inspired_theme()

        # todo(rteqs): figure out why we can't just catch KeyboardInterrupt without crashing the kernel
        signal.signal(
            signal.SIGUSR1,
            lambda signum, frame: cell_interrupt()
            if self.active_cell is not None
            and self.cell_status[self.active_cell] == "running"
            else None,
        )

    def debug_state(self) -> dict[str, object]:
        return {
            "cell_seq": self.cell_seq,
            "cell_rnodes": {k: v.debug_state() for k, v in self.cell_rnodes.items()},
            "k_globals": {
                "touched": list(self.k_globals.touched),
                "removed": list(self.k_globals.removed),
                "dataframes": list(self.k_globals.dataframes.sample()),
            },
            "cell_status": self.cell_status,
            "active_cell": self.active_cell,
            "widget_signals": {k: repr(v) for k, v in self.widget_signals.items()},
            "nodes_with_widgets": {
                str(k): v.debug_state() for k, v in self.nodes_with_widgets.items()
            },
            "cell_output_selections": self.cell_output_selections,
            "viewer_cell_selections": self.viewer_cell_selections,
            "plot_data_selections": self.plot_data_selections,
            "cell_pagination_settings": {
                cell_id: {k: asdict(settings) for k, settings in d.items()}
                for cell_id, d in self.cell_pagination_settings.items()
            },
            "viewer_pagination_settings": {
                viewer_id: {k: asdict(settings) for k, settings in d.items()}
                for viewer_id, d in self.viewer_pagination_settings.items()
            },
            "ldata_dataframes": list(self.ldata_dataframes.keys()),
            "registry_dataframes": list(self.registry_dataframes.keys()),
            "url_dataframes": list(self.url_dataframes.keys()),
            "plot_configs": self.plot_configs,
        }

    async def set_active_cell(self, cell_id: str) -> None:
        self.active_cell = cell_id

        # todo(maximsmol): huge hack to make sure runtime has time to read
        # our std streams before we ask for active cell to be switched
        #
        # we need to instead overwrite stdout/stderr with a thing that appends
        # headers indicating the active cell
        sys.stdout.flush()
        sys.stderr.flush()
        await asyncio.sleep(0.1)

        self.cell_seq += 1
        await self.send(
            {"type": "start_cell", "cell_id": cell_id, "run_sequencer": self.cell_seq}
        )

    async def send_global_updates(self) -> None:
        async with asyncio.TaskGroup() as tg:
            for cell_id, key in self.cell_output_selections.items():
                if key not in self.k_globals.touched:
                    continue

                tg.create_task(self.send_output_value(key, cell_id=cell_id))

            # note: plots with filter tables have two layers of indirection.
            # key -> controlling viewer -> viewer
            touched_viewers = {
                f"df_{viewer_id}"
                for viewer_id, sources in self.viewer_pagination_settings.items()
                if any(key in self.k_globals.touched for key in sources)
            }
            touched_viewers |= {
                f"df_{viewer_id}"
                for viewer_id, sources in self.viewer_pagination_settings.items()
                if any(key in touched_viewers for key in sources)
            }

            for viewer_id, (key, key_type) in self.viewer_cell_selections.items():
                if key not in self.k_globals.touched and key not in touched_viewers:
                    continue

                if key_type == "key":
                    tg.create_task(self.send_output_value(key, viewer_id=viewer_id))
                elif key_type == "ldata_node_id":
                    tg.create_task(
                        self.send_output_value(ldata_node_id=key, viewer_id=viewer_id)
                    )
                elif key_type == "registry_table_id":
                    tg.create_task(
                        self.send_output_value(
                            registry_table_id=key, viewer_id=viewer_id
                        )
                    )
                elif key_type == "url":
                    tg.create_task(self.send_output_value(url=key, viewer_id=viewer_id))

            for plot_id, key in self.plot_data_selections.items():
                if key not in self.k_globals.touched and key not in touched_viewers:
                    continue

                tg.create_task(self.send_plot_data(plot_id, key))

            tg.create_task(self.send_globals_summary())

    async def on_tick_finished(self, updated_signals: dict[int, Signal]) -> None:
        # todo(maximsmol): this can be optimizied
        # 1. we can just update nodes that actually re-ran last tick instead of everything
        # 2. we can pre-compute nww_by_cell

        nww_by_cell = {
            cell_id: [
                x for x in self.nodes_with_widgets.values() if x.cell_id == cell_id
            ]
            for cell_id in self.cell_rnodes
        }

        updated_widgets: set[str] = set()

        unused_signals: set[str] = set(self.widget_signals.keys())
        for cell_id in self.cell_rnodes:
            res: dict[str, WidgetState] = {}

            for n in nww_by_cell[cell_id]:
                path = n.name_path()
                for k, v in n.widget_states.items():
                    abs_k = f"{path}/{k}"
                    res[abs_k] = {**v}

                    if abs_k not in self.widget_signals:
                        continue
                    unused_signals.remove(abs_k)

                    sig = self.widget_signals[abs_k]
                    if id(sig) in updated_signals:
                        updated_widgets.add(abs_k)

                    val = sig.sample()
                    if val is Nothing.x:
                        continue

                    res[abs_k]["value"] = val
            if self.cell_status[cell_id] == "error":
                # skip errored cells to avoid clobbering widget state
                # must be here so that we update unused_signals properly
                # todo(maximsmol): optimize
                continue

            await self.send(
                {
                    "type": "cell_widgets",
                    "cell_id": cell_id,
                    "widget_state": res,
                    "updated_widgets": list(updated_widgets),
                }
            )

        # fixme(rteqs): cleanup signals in some other way. the below does not work because widget signals
        # are restored on `init` but there are no corresponding `rnodes`
        # for x in unused_signals:
        #     del self.widget_signals[x]

    def get_widget_value(self, key: str) -> Signal[Any]:
        assert ctx.cur_comp is not None

        key = f"{ctx.cur_comp.name_path()}/{key}"

        if key not in self.widget_signals:
            self.widget_signals[key] = Signal(Nothing.x, name=key)

        return self.widget_signals[key]

    def emit_widget(self, key: str, data: WidgetState) -> None:
        assert ctx.cur_comp is not None

        ctx.cur_comp.widget_states[key] = data
        self.nodes_with_widgets[id(ctx.cur_comp)] = ctx.cur_comp

    def submit_widget_state(self) -> None:
        for s in ctx.updated_signals.values():
            s._apply_updates()

        self.conn.call_fut(self.on_tick_finished(ctx.signals_update_from_code)).result()

    def on_dispose(self, node: Node) -> None:
        if id(node) not in self.nodes_with_widgets:
            return

        del self.nodes_with_widgets[id(node)]

    def _cell_outputs(self) -> CategorizedCellOutputs:
        res = CategorizedCellOutputs()
        for x in self.k_globals.available:
            res.all.append(x)

            val = self.k_globals[x]
            if hasattr(val, "iloc"):
                res.dfs.append(x)

            if isinstance(val, BaseFigure):
                res.figures.append(x)

            elif isinstance(val, Figure) or (
                hasattr(val, "figure") and isinstance(val.figure, Figure)
            ):
                res.static_figures.append(x)

        res.all.sort()
        res.dfs.sort()
        res.figures.sort()
        res.static_figures.sort()

        return res

    async def exec(self, *, cell_id: str, code: str) -> None:
        filename = f"<cell {cell_id}>"

        try:
            assert ctx.cur_comp is None
            assert not ctx.in_tx

            self.cell_status[cell_id] = "running"

            comp = self.cell_rnodes.get(cell_id)
            if comp is not None:
                comp.dispose()
                del self.cell_rnodes[cell_id]

            # https://stackoverflow.com/questions/33908794/get-value-of-last-expression-in-exec-call
            parsed = compile(
                source=code,
                filename=filename,
                mode="exec",
                flags=ast.PyCF_ONLY_AST | ast.PyCF_ALLOW_TOP_LEVEL_AWAIT,
            )

            stmts = list(ast.iter_child_nodes(parsed))
            if len(stmts) == 0:
                self.cell_status[cell_id] = "ok"
                self.k_globals.clear()
                await self.send_cell_result(cell_id)
                return

            async def x() -> None:
                try:
                    assert ctx.cur_comp is not None

                    self.cell_rnodes[cell_id] = ctx.cur_comp
                    self.k_globals.clear()

                    try:
                        res = eval(  # noqa: S307
                            compile(
                                parsed,
                                filename=filename,
                                mode="exec",
                                flags=ast.PyCF_ALLOW_TOP_LEVEL_AWAIT,
                            ),
                            self.k_globals,
                        )
                        if asyncio.iscoroutine(res):
                            res = await res
                    except ExitException:
                        ...

                    self.cell_status[cell_id] = "ok"
                    await self.send_cell_result(cell_id)

                except Exception:
                    self.cell_status[cell_id] = "error"
                    await self.send_cell_result(cell_id)

            x.__name__ = filename

            await self.set_active_cell(cell_id)
            await ctx.run(x, _cell_id=cell_id)

        except Exception:
            self.cell_status[cell_id] = "error"
            await self.send_cell_result(cell_id)

    async def send_cell_result(self, cell_id: str) -> None:
        outputs = sorted(self.k_globals.available)
        outputs = self._cell_outputs()

        msg = {
            "type": "cell_result",
            "cell_id": cell_id,
            "outputs": outputs.all,
            "dataframe_outputs": outputs.dfs,
            "figure_outputs": outputs.figures,
            "static_figure_outputs": outputs.static_figures,
        }
        if sys.exception() is not None:
            msg["exception"] = format_exc()

        await self.send(msg)
        await self.send_global_updates()

    async def send(self, msg: object) -> None:
        # print("[kernel] >", msg)
        await self.conn.send(msg)

    async def send_plot_data(
        self, plot_id: str, key: str, config: PlotConfig | None = None
    ) -> None:
        if config is not None and config == self.plot_configs.get(plot_id):
            return

        self.plot_configs[plot_id] = config

        res = self.k_globals[key]

        if isinstance(res, BaseFigure):
            await self.send(
                {
                    "type": "plot_data",
                    "plot_id": plot_id,
                    "key": key,
                    # todo(maximsmol): get rid of the json reload
                    "plotly_json": serialize_plotly_figure(res),
                }
            )
            return

        if not hasattr(res, "iloc"):
            await self.send(
                {
                    "type": "plot_data",
                    "plot_id": plot_id,
                    "key": key,
                    "error": "not a dataframe",
                }
            )
            return

        if hasattr(res, "compute"):
            # Dask support
            res = res.compute()

        if isinstance(res, Series):
            res = pd.DataFrame(res)

        if isinstance(res, DataFrame):
            msg = {
                "type": "plot_data",
                "plot_id": plot_id,
                "key": key,
                "dataframe_json": {"schema": build_table_schema(res, version=False)},
            }

            df_size_mb = res.memory_usage(index=True, deep=True).sum() / 10**6

            if df_size_mb <= 10:
                # todo(maximsmol): get rid of the json reload
                msg["dataframe_json"]["data"] = orjson.loads(
                    res.to_json(orient="split", date_format="iso")
                )

            elif self.duckdb is not None and config is not None:
                subsampled_data = await downsample_df(
                    self.duckdb,
                    key,
                    res,
                    config,
                    self.k_globals.item_write_counter[key],
                )
                msg["dataframe_json"]["subsampled_data"] = subsampled_data

                for trace in config.get("traces", []):
                    if trace["type"] != "scattergl" and trace["type"] != "scatter":
                        msg["dataframe_json"]["data"] = orjson.loads(
                            res.to_json(orient="split", date_format="iso")
                        )
                        break

            await self.send(msg)

    def lookup_pagination_settings(
        self,
        *,
        cell_id: str | None = None,
        viewer_id: str | None = None,
        source_key: str,
    ) -> PaginationSettings:
        if cell_id is not None:
            return self.cell_pagination_settings[cell_id][source_key]
        if viewer_id is not None:
            return self.viewer_pagination_settings[viewer_id][source_key]

        raise ValueError("Both cell_id and viewer_id are null")

    async def send_output_value(
        self,
        key: str | None = None,
        ldata_node_id: str | None = None,
        registry_table_id: str | None = None,
        url: str | None = None,
        *,
        cell_id: str | None = None,
        viewer_id: str | None = None,
    ) -> None:
        res: Any
        key_fields: dict[str, str] | None = None

        if key is not None:
            res = self.k_globals[key]
            key_fields = {"key": key}
        elif ldata_node_id is not None:
            res = self.ldata_dataframes[ldata_node_id]
            key_fields = {"ldata_node_id": ldata_node_id}
        elif registry_table_id is not None:
            res = self.registry_dataframes[registry_table_id]
            key_fields = {"registry_table_id": registry_table_id}
        elif url is not None:
            res = self.url_dataframes[url]
            key_fields = {"url": url}
        else:
            raise RuntimeError(
                "Requested value without key, ldata_node_id, registry_table_id or url"
            )

        id_fields: dict[str, str] | None = None
        if cell_id is not None:
            id_fields = {"cell_id": cell_id}

        if viewer_id is not None:
            id_fields = {"viewer_id": viewer_id}

        assert id_fields is not None
        assert key_fields is not None
        assert len(key_fields) == 1

        if isinstance(res, BaseFigure):
            await self.send(
                {
                    "type": "output_value",
                    **(id_fields),
                    **(key_fields),
                    "plotly_json": orjson.dumps(
                        serialize_plotly_figure(res), option=orjson.OPT_SERIALIZE_NUMPY
                    ).decode(),
                }
            )
            return

        if hasattr(res, "iloc"):
            if isinstance(res, Series):
                res = pd.DataFrame(res)

            pagination_settings = self.lookup_pagination_settings(
                cell_id=cell_id,
                viewer_id=viewer_id,
                source_key=next(iter(key_fields.values())),
            )
            sort_settings = pagination_settings.sort_settings
            row_filters = pagination_settings.row_filters

            res = filter_and_sort(df=res, pagination_settings=pagination_settings)

            if viewer_id is not None:
                async with ctx.transaction:
                    filtered_dataframe_name = f"df_{viewer_id}"
                    self.k_globals[filtered_dataframe_name] = res

            page_size = pagination_settings.page_size
            page_idx = pagination_settings.page_idx

            num_pages = math.ceil(len(res) / page_size)
            h = page_size

            data = paginate(df=res, pagination_settings=pagination_settings)

            # we don't want to repeat the column names
            # so we use the "split" format
            # but we still want the schema
            await self.send(
                {
                    "type": "output_value",
                    **(id_fields),
                    **(key_fields),
                    "dataframe_json": {
                        "schema": build_table_schema(data, version=False),
                        "data": data.to_dict(orient="split"),
                        # todo(maximsmol): this seems useless?
                        "num_pages": num_pages,
                        "page_idx": page_idx,
                        "page_size": h,
                        **(
                            {"sort_settings": sort_settings}
                            if sort_settings is not None
                            else {}
                        ),
                        **(
                            {"row_filters": row_filters}
                            if row_filters is not None
                            else {}
                        ),
                    },
                }
            )
            return

        if hasattr(res, "__dataframe__"):
            # dataframe interchange object
            # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.__dataframe__.html
            # there is a lot to implement here so it's not really supported right now
            # e.g. we need DLPack to properly do this:
            # https://dmlc.github.io/dlpack/latest/python_spec.html
            # tho pandas itself does not support DLPack yet
            # reference implementation:
            # https://github.com/pandas-dev/pandas/blob/v2.2.2/pandas/core/interchange/from_dataframe.py#L33-L73

            data = res.__dataframe__()

            await self.send(
                {
                    "type": "output_value",
                    **(id_fields),
                    **(key_fields),
                    "__dataframe__": {
                        "num_columns": data.num_columns(),
                        "num_rows": data.num_rows(),
                    },
                }
            )
            return

        if isinstance(res, Figure) or (
            hasattr(res, "figure") and isinstance(res.figure, Figure)
        ):
            await self.send(
                {"type": "output_value", **(id_fields), **(key_fields), "webp": res}
            )
            return

        data = pprint.pformat(res)
        await self.send(
            {
                "type": "output_value",
                **(id_fields),
                **(key_fields),
                "string": data[:10000],
            }
        )

    async def send_globals_summary(self) -> None:
        summary = {}
        for key, value in self.k_globals.items():
            if isinstance(value, Signal):
                value = value.sample()

            if isinstance(value, pd.DataFrame):
                summary[key] = {
                    "type": "DataFrame",
                    "columns": list(value.columns)[:50],
                    "dtypes": {
                        str(k): str(v) for k, v in list(value.dtypes.items())[:50]
                    },
                    "column_preview_truncated": len(value.columns) > 50,
                    "shape": value.shape,
                }
            elif isinstance(value, pd.Series):
                summary[key] = {
                    "type": "Series",
                    "dtype": str(value.dtype),
                    "shape": value.shape,
                }
            else:
                summary[key] = {"type": type(value).__name__}

        await self.send({"type": "globals_summary", "summary": summary})

    async def accept(self) -> None:
        # print("[kernel] accept")
        msg = await self.conn.recv()
        # print("[kernel] <", msg)

        if msg["type"] == "init":
            self.cell_output_selections = msg["cell_output_selections"]
            self.plot_data_selections = msg["plot_data_selections"]
            self.plot_configs = msg["plot_configs"]

            viewer_cell_data = msg["viewer_cell_data"]
            for cell_id, data in viewer_cell_data.items():
                self.viewer_cell_selections[cell_id] = (
                    data["source"][0],
                    data["source"][1],
                )

                for k, v in data["pagination_settings"].items():
                    setattr(
                        self.viewer_pagination_settings[cell_id][data["source"][0]],
                        k,
                        v,
                    )

            for state_raw in msg["widget_states"].values():
                try:
                    state: dict[str, WidgetState] = orjson.loads(state_raw)
                except orjson.JSONDecodeError:
                    continue

                for k, v in state.items():
                    if "value" not in v:
                        continue

                    self.widget_signals[k] = Signal(v["value"])

            async with ctx.transaction:
                for viewer_id, (
                    data_id,
                    key_type,
                ) in self.viewer_cell_selections.items():
                    key = f"df_{viewer_id}"

                    df: DataFrame | None = None
                    if key_type == "ldata_node_id":
                        path_str = f"latch://{data_id}.node"
                        presigned_url = await get_presigned_url(path_str)

                        if data_id not in self.ldata_dataframes:
                            self.ldata_dataframes[data_id] = pd.read_csv(
                                presigned_url, sep=None, engine="python"
                            )

                        df = self.ldata_dataframes[data_id]

                    elif key_type == "registry_table_id":
                        table = Table(id=data_id)
                        table_df = table.get_dataframe()

                        if data_id not in self.registry_dataframes:
                            self.registry_dataframes[data_id] = table_df

                        df = self.registry_dataframes[data_id]

                    elif key_type == "url":
                        if data_id not in self.url_dataframes:
                            self.url_dataframes[data_id] = pd.read_csv(
                                data_id, sep=None, engine="python"
                            )

                        df = self.url_dataframes[data_id]

                    if df is not None:
                        pagination_settings = self.lookup_pagination_settings(
                            viewer_id=viewer_id, source_key=data_id
                        )
                        self.k_globals[key] = filter_and_sort(
                            df=df, pagination_settings=pagination_settings
                        )

            return

        if msg["type"] == "debug_state":
            await self.send({"type": "debug_state", "data": self.debug_state()})
            return

        if msg["type"] == "run_cell":
            await self.exec(cell_id=msg["cell_id"], code=msg["code"])
            return

        if msg["type"] == "dispose_cell":
            cell_id = msg["cell_id"]

            node = self.cell_rnodes.get(cell_id)
            if node is not None:
                node.dispose()
                del self.cell_rnodes[cell_id]
                del self.cell_status[cell_id]

            return

        if msg["type"] == "get_plot_data":
            self.plot_data_selections[msg["plot_id"]] = msg["key"]

            await self.send_plot_data(msg["plot_id"], msg["key"], msg.get("config"))

        if msg["type"] == "get_global":
            pagination_settings = {
                setting_name: msg.get(setting_name, default)
                for setting_name, default in {
                    "page_size": 25,
                    "page_idx": 0,
                    "row_filters": [],
                    "sort_settings": None,
                    "selections": None,
                }.items()
            }

            if "cell_id" in msg:
                self.cell_output_selections[msg["cell_id"]] = msg["key"]
                for k, v in pagination_settings.items():
                    setattr(
                        self.cell_pagination_settings[msg["cell_id"]][msg["key"]], k, v
                    )
                await self.send_output_value(msg["key"], cell_id=msg["cell_id"])

            if "viewer_id" in msg:
                viewer_id = msg["viewer_id"]

                key_fields: dict[KeyType, str] | None = None

                if "key" in msg:
                    key_fields = {"key": msg["key"]}

                elif "ldata_node_id" in msg:
                    ldata_node_id = msg["ldata_node_id"]
                    path_str = f"latch://{ldata_node_id}.node"

                    if ldata_node_id not in self.ldata_dataframes:
                        presigned_url = await get_presigned_url(path_str)
                        self.ldata_dataframes[ldata_node_id] = pd.read_csv(
                            presigned_url, sep=None, engine="python"
                        )

                    key_fields = {"ldata_node_id": ldata_node_id}

                elif "registry_table_id" in msg:
                    registry_table_id = msg["registry_table_id"]
                    table = Table(id=registry_table_id)
                    table_df = table.get_dataframe()

                    if registry_table_id not in self.registry_dataframes:
                        self.registry_dataframes[registry_table_id] = table_df

                    key_fields = {"registry_table_id": registry_table_id}

                elif "url" in msg:
                    url = msg["url"]
                    if url not in self.url_dataframes:
                        self.url_dataframes[url] = pd.read_csv(
                            url, sep=None, engine="python"
                        )

                    key_fields = {"url": url}

                else:
                    raise RuntimeError(
                        "Requested value without key, ldata_node_id, registry_table_id or url"
                    )

                assert key_fields is not None
                assert len(key_fields) == 1

                key_type, key = next(iter(key_fields.items()))
                self.viewer_cell_selections[viewer_id] = (key, key_type)

                for k, v in pagination_settings.items():
                    setattr(
                        self.viewer_pagination_settings[msg["viewer_id"]][key], k, v
                    )

                await self.send_output_value(
                    **(key_fields), viewer_id=msg.get("viewer_id")
                )

                await self.send_global_updates()

        if msg["type"] == "set_widget_value":
            for w_key, payload in msg["data"].items():
                try:
                    if w_key not in self.widget_signals:
                        continue

                    async with ctx.transaction:
                        self.widget_signals[w_key](
                            orjson.loads(payload), _ui_update=True
                        )
                except Exception:
                    traceback.print_exc()
                    continue

            return

        if msg["type"] == "globals_summary":
            await self.send_globals_summary()
            return


loop: asyncio.AbstractEventLoop | None = None
shutdown_requested = False


def sigterm_handler(signum: int, frame: FrameType | None) -> None:
    for t in asyncio.all_tasks(loop):
        t.cancel("SIGTERM")

    # try to exit using KeyboardInterrupt
    signal.raise_signal(signal.SIGINT)


async def main() -> None:
    global loop
    loop = asyncio.get_running_loop()

    signal.signal(signal.SIGTERM, sigterm_handler)

    sock = socket.socket(family=socket.AF_UNIX, fileno=int(sys.argv[-1]))
    sock.setblocking(False)

    # note(maximsmol): DO NOT LET THESE GET GARBAGE COLLECTED
    # the interpreter will crash if that happens apparently
    old_stdout = sys.stdout
    old_stderr = sys.stderr

    socket_io_thread = SocketIoThread(socket=sock)
    socket_io_thread.start()
    try:
        socket_io_thread.initialized.wait()

        k = Kernel(conn=socket_io_thread)
        _inject.kernel = k

        sys.stdout = text_socket_writer(
            SocketWriter(conn=k.conn, kernel=k, name="stdout", loop=loop)
        )
        sys.stderr = text_socket_writer(
            SocketWriter(conn=k.conn, kernel=k, name="stderr", loop=loop)
        )

        await k.send({"type": "ready"})

        while not shutdown_requested:
            try:
                await k.accept()
            except Exception:
                traceback.print_exc()
                continue

        print("Kernel shutting down...")
    finally:
        socket_io_thread.shutdown.set()
        socket_io_thread.join()

        sys.stdout = old_stdout
        sys.stderr = old_stderr


if __name__ == "__main__":
    import multiprocessing

    multiprocessing.set_start_method("forkserver")

    asyncio.run(main())
