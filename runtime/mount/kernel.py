if __name__ == "__main__":
    # fixme(maximsmol): anndata calls multiprocessing.get_context
    # so we have to do this BEFORE we import it indirectly
    import multiprocessing

    multiprocessing.set_start_method("forkserver")

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
from copy import copy, deepcopy
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
from latch.ldata.path import LPath
from latch.registry.table import Table
from latch_cli.utils import urljoins
from lplots import _inject
from lplots.h5.process_message import handle_h5_widget_message
from lplots.h5.utils import auto_install
from lplots.h5.utils.persistence import load_anndata, serialize_anndata
from lplots.persistence import (
    SerializedNode,
    SerializedSignal,
    safe_serialize_obj,
    safe_unserialize_obj,
    small_repr,
    unable_to_unserialize_symbol,
)
from lplots.reactive import Node, Signal, ctx, live_nodes, live_signals
from lplots.themes import graphpad_inspired_theme
from lplots.utils.nothing import Nothing
from lplots.widgets._emit import WidgetState
from lplots.widgets.widget import BaseWidget, load_widget_helper
from matplotlib.figure import Figure
from pandas import DataFrame, Index, MultiIndex, Series
from pandas.io.json._table_schema import build_table_schema
from plotly.basedatatypes import BaseFigure
from plotly_utils.precalc_box import precalc_box
from plotly_utils.precalc_violin import precalc_violin
from socketio_thread import SocketIoThread
from stdio_over_socket import SocketWriter, text_socket_writer

ad = auto_install.ad

sys.path.append(str(Path(__file__).parent.absolute()))
from subsample import downsample_df, initialize_duckdb
from utils import KernelSnapshotStatus, PlotConfig, get_presigned_url

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


class TracedDict(dict[str, Signal[object] | object]):
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

    def _direct_set(self, __key: str, __value: object) -> None:
        return super().__setitem__(__key, __value)

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


class ExitException(Exception):
    ...


KeyType = Literal["key", "ldata_node_id", "registry_table_id", "url"]


def cell_exit(code: int = 0) -> None:
    raise ExitException


def cell_interrupt(code: int = 0) -> None:
    raise KeyboardInterrupt


leading_digits_and_dash = re.compile(r"^\d+-")


def remove_leading_digits_and_dash(col_raw: str) -> str:
    return re.sub(leading_digits_and_dash, "", col_raw)


multi_index_col_name = re.compile(r"^level_\d+$")


def is_multi_index_col(col: str) -> bool:
    return re.match(multi_index_col_name, col) is not None


K = TypeVar("K")
V = TypeVar("V")
U = TypeVar("U")


def filter_items(source: dict[K, V], lookup: dict[V, U]) -> dict[K, U]:
    return {k: value for k, v in source.items() if (value := lookup.get(v)) is not None}


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

    mask: Series[bool] | NDArray[np.bool_] | None = None

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
) -> "DataFrame | Series[bool] | NDArray[np.bool_]":
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
        data = df.iloc[:, x : x + w].head(y + h, npartitions=-1)[y:]
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

def _split_violin_groups(trace: dict[str, Any]) -> list[dict[str, Any]] | None:
    orientation = trace.get("orientation", "v")
    data_axis = "y" if orientation == "v" else "x"
    index_axis = "x" if orientation == "v" else "y"
    if index_axis not in trace or data_axis not in trace:
        return None

    # NOTE(tim): handle binary plotly typed array that's in the trace
    def _decode_plotly_typed(val: Any) -> Any:
        if isinstance(val, dict) and "bdata" in val and "dtype" in val:
            import base64 as _b64
            buf = _b64.b64decode(val["bdata"])
            return np.frombuffer(buf, dtype=np.dtype(val["dtype"]))
        return val

    idx_arr = np.asarray(trace.get(index_axis))
    vals_arr = np.asarray(_decode_plotly_typed(trace.get(data_axis)))

    print('idx_arr', idx_arr, 'vals_arr', vals_arr)
    print('getattr(idx_arr, "ndim", 1)', getattr(idx_arr, "ndim", 1), 'getattr(vals_arr, "ndim", 1)', getattr(vals_arr, "ndim", 1))
    # TODO(tim): consider deleting this err handling
    # make sure input is 1D, not empty, and has the same length
    if getattr(idx_arr, "ndim", 1) != 1 or getattr(vals_arr, "ndim", 1) != 1:
        return None
    n = int(vals_arr.shape[0])
    if n == 0 or int(idx_arr.shape[0]) != n:
        return None

    # group by category using NumPy with sorted category order
    categories, cat_idx = np.unique(idx_arr, return_inverse=True)
    if len(categories) <= 1:
        return None
    order = categories.tolist()
    cat_idx = cat_idx.astype(np.int64, copy=False)
    order_idx = np.argsort(cat_idx, kind="stable")
    cat_idx_sorted = cat_idx[order_idx]
    # find where to split the categories
    split_idx = np.flatnonzero(np.diff(cat_idx_sorted)) + 1
    groups = np.split(order_idx, split_idx)
    print('groups_per_category', [g.size for g in groups], 'order', order)
    group_traces: list[dict[str, Any]] = []
    for label, idxs in zip(order, groups):
        child = deepcopy(trace)
        child[data_axis] = vals_arr[idxs]
        # child[index_axis] = [label] * len(child[data_axis])
        # clear any indexes that will cause trouble in precalc_violin
        child.pop(index_axis, None)
        child.pop(f"{index_axis}0", None)
        child.pop(f"d{index_axis}", None)
        child["name"] = str(label)
        print(
            "[split_violin] child:",
            "label=", label,
            "name=", child.get("name"),
            "n=", int(idxs.size),
        )
        group_traces.append(child)
    return group_traces


def serialize_plotly_figure(x: BaseFigure) -> object:
    res = x.to_dict()
    orig_n = len(res["data"]) if isinstance(res.get("data"), list) else 0

    data_out: list[dict[str, Any]] = []
    for trace in res["data"]:
        try:
            if trace["type"] == "box":
                precalc_box(trace)
                data_out.append(trace)
            elif trace["type"] == "violin":
                # NOTE(tim): if the trace has multiple violins,
                # seperate them into seperate traces so we can precompute
                # them separately
                group_traces = _split_violin_groups(trace)
                if group_traces is not None:
                    try:
                        print(
                            f"[plots-faas] violin fan-out: groups={len(group_traces)} name={trace.get('name')} orient={trace.get('orientation','v')}",
                            flush=True,
                        )
                    except Exception:
                        pass
                    for group_trace in group_traces:
                        try:
                            print("group trace before precalc", group_trace)

                            

                            precalc_violin(group_trace)
                            print("group trace after precalc", group_trace)

                            # Position the split violin at its category via a
                            # single-element position array (e.g., x=[label]) and
                            # clear any conflicting anchors/offsets so the label
                            # aligns with center.

                            orientation = trace.get("orientation", "v")
                            data_axis = "y" if orientation == "v" else "x"
                            index_axis = "x" if orientation == "v" else "y"
                            label     = str(group_trace.get("name", ""))
                            # group_trace[index_axis] = [label] * 10
                            # NOTE: this col is not being used, but plotly has 
                            # weird issues if it's not provided, so set it to 0
                            group_trace[data_axis] = []
                            # group_trace[data_axis] = [label] * 10
                            # group_trace[pos_axis] = [label] * 10
                            # group_trace[data_axis] = [label]
                            # group_trace[pos_axis] = [label] * len(group_trace[data_axis])
                            # group_trace[data_axis] = [label]
                            # print(f"group_trace[data_axis] {group_trace[data_axis]}")
                            # group_trace[pos_axis] = [label] * len(group_trace[data_axis])
                            # group_trace.pop(f"{pos_axis}0", None)
                            # group_trace.pop(f"d{pos_axis}", None)
                            # group_trace.pop(f"{data_axis}0", None)
                            # group_trace.pop(f"d{data_axis}", None)
                            # group_trace.pop("offsetgroup", None)
                            print(f"group trace after post compute {group_trace}")

                        except Exception:
                            traceback.print_exc()
                        data_out.append(group_trace)
                else:
                    precalc_violin(trace)
                    data_out.append(trace)
            elif trace["type"] == "scatter":
                trace["type"] = "scattergl"
                data_out.append(trace)
            else:
                data_out.append(trace)
        except Exception:
            traceback.print_exc()
            data_out.append(trace)

    res["data"] = data_out
    try:
        print(
            f"[plots-faas] serialize_plotly_figure: traces in={orig_n}, out={len(res['data'])}",
            flush=True,
        )
    except Exception:
        pass

    modules = {
        "sage_all": pio_json.get_module("sage.all", should_load=False),
        "np": pio_json.get_module("numpy", should_load=False),
        "pd": pio_json.get_module("pandas", should_load=False),
        "image": pio_json.get_module("PIL.Image", should_load=False),
    }

    # note(maximsmol): plotly itself does a bunch of escaping to avoid XSS
    # when embedding directly into HTML. we never do that so we don't care
    return pio_json.clean_to_json_compatible(res, modules=modules)


snapshot_dir = Path.home() / ".cache" / "plots-faas"
snapshot_f_name = "snapshot.json"


class SerializedGlobal(TypedDict):
    value: str
    error_msg: str | None


class RestoredGlobalInfo(TypedDict):
    value: object
    msg: str


snapshot_chunk_bytes = 64 * 2 ** 10
snapshot_progress_interval_bytes = 10 * (2 ** 20)


@dataclass(kw_only=True)
class Kernel:
    conn: SocketIoThread

    cell_seq: int = 0
    cell_rnodes: dict[str, Node] = field(default_factory=dict)
    k_globals: TracedDict = field(init=False)
    cell_status: dict[str, str] = field(default_factory=dict)
    snapshot_status: KernelSnapshotStatus = "done"

    active_cell: str | None = None

    widget_signals: dict[str, Signal[Any]] = field(default_factory=dict)
    nodes_with_widgets: dict[str, Node] = field(default_factory=dict)

    cells_with_pending_widget_updates: set[str] = field(default_factory=set)

    cell_output_selections: dict[str, str] = field(default_factory=dict)
    viewer_cell_selections: dict[str, tuple[str, KeyType]] = field(default_factory=dict)
    plot_data_selections: dict[str, str] = field(default_factory=dict)

    ldata_dataframes: dict[str, DataFrame] = field(default_factory=dict)
    registry_dataframes: dict[str, DataFrame] = field(default_factory=dict)
    url_dataframes: dict[str, DataFrame] = field(default_factory=dict)

    ann_data_objects: dict[str, ad.AnnData] = field(default_factory=dict)

    cell_pagination_settings: defaultdict[str, defaultdict[str, PaginationSettings]] = (
        field(default_factory=pagination_settings_dict_factory)
    )
    viewer_pagination_settings: defaultdict[
        str, defaultdict[str, PaginationSettings]
    ] = field(default_factory=pagination_settings_dict_factory)

    plot_configs: dict[str, PlotConfig | None] = field(default_factory=dict)
    duckdb: DuckDBPyConnection = field(default=initialize_duckdb())

    session_snapshot_mode: bool = False

    restored_nodes: dict[str, Node] = field(default_factory=dict)
    restored_signals: dict[str, Signal[object]] = field(default_factory=dict)
    restored_globals: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.k_globals = TracedDict(self.duckdb)
        self.k_globals["exit"] = cell_exit
        self.k_globals.clear()
        pio.templates["graphpad_inspired_theme"] = graphpad_inspired_theme()

        signal.signal(
            signal.SIGINT,
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
            "widget_signals": {
                k: v.serialize(short_val=True) for k, v in self.widget_signals.items()
            },
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
            "restored_nodes": {
                k: v.serialize() for k, v in self.restored_nodes.items()
            },
            "restored_signals": {
                k: v.serialize(short_val=True) for k, v in self.restored_signals.items()
            },
            "live_signals": {k: v.serialize() for k, v in live_signals.items()},
            "live_nodes": {k: v.serialize() for k, v in live_nodes.items()},
            "restored_globals": self.restored_globals,
            "session_snapshot_mode": self.session_snapshot_mode,
        }

    # note(maximsmol): called by the reactive context
    async def set_active_cell(self, cell_id: str | None) -> None:
        # todo(maximsmol): I still believe this is correct
        # but we need to deal with the frontend clearing logs in weird ways
        # e.g. a button should maybe clear the logs if it triggered its own cell
        # but a subnode should probably not clear logs?
        #
        # right now we rely on start_cell to clear logs each time anything in a cell runs

        # if self.active_cell == cell_id:
        #     return

        # note(maximsmol): stdio_over_socket will fetch the active cell id
        # synchronously in `.write` (which will be called by the buffered writers' `.flush`)
        sys.stdout.flush()
        sys.stderr.flush()

        if cell_id is not None:
            self.active_cell = cell_id

        self.cell_seq += 1
        await self.send(
            {
                "type": "start_cell",
                "cell_id": cell_id,
                "run_sequencer": self.cell_seq,
            }
        )

    async def send_global_updates(self) -> None:
        async with asyncio.TaskGroup() as tg:
            for cell_id, key in self.cell_output_selections.items():
                if key not in self.k_globals.available:
                    continue

                tg.create_task(self.send_output_value(key, cell_id=cell_id))

            # note: plots with filter tables have two layers of indirection.
            # key -> controlling viewer -> viewer
            touched_viewers = {
                f"df_{viewer_id}"
                for viewer_id, sources in self.viewer_pagination_settings.items()
                if any(key in self.k_globals.available for key in sources)
            }
            touched_viewers |= {
                f"df_{viewer_id}"
                for viewer_id, sources in self.viewer_pagination_settings.items()
                if any(key in touched_viewers for key in sources)
            }

            for viewer_id, (key, key_type) in self.viewer_cell_selections.items():
                if key not in self.k_globals.available and key not in touched_viewers:
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
                if key not in self.k_globals.available and key not in touched_viewers:
                    continue

                tg.create_task(self.send_plot_data(plot_id, key))

            tg.create_task(self.send_globals_summary())

    async def on_tick_finished(
        self, updated_signals: dict[int, Signal[object]], clear_status: bool = True
    ) -> None:
        # todo(maximsmol): this can be optimizied
        # 1. we can just update nodes that actually re-ran last tick instead of everything
        # 2. we can pre-compute nww_by_cell
        # 3. we can batch together the update messages so the frontend does not
        # need to spam GQL requests

        nww_by_cell = {
            cell_id: [
                x for x in self.nodes_with_widgets.values() if x.cell_id == cell_id
            ]
            for cell_id in self.cell_rnodes
        }

        for cell_id in self.cell_rnodes:
            # note(kenny): if kernel loaded from snapshot, not every cell rnode
            # will have a status
            if self.cell_status.get(cell_id) == "error":
                # skip errored cells to avoid clobbering widget state
                continue

            updated_widgets: set[str] = set()
            res: dict[str, WidgetState[str, object]] = {}
            for n in nww_by_cell[cell_id]:
                path = n.name_path()
                for k, v in n.widget_states.items():
                    abs_k = f"{path}/{k}"
                    res[abs_k] = copy(v)

                    if abs_k not in self.widget_signals:
                        continue

                    sig = self.widget_signals[abs_k]
                    if sig.id in updated_signals:
                        updated_widgets.add(abs_k)

                    val = sig.sample()
                    if val is Nothing.x:
                        continue

                    res[abs_k]["value"] = val

            if (
                len(updated_widgets) == 0
                and cell_id not in self.cells_with_pending_widget_updates
            ):
                continue

            await self.send(
                {
                    "type": "cell_widgets",
                    "cell_id": cell_id,
                    "widget_state": res,
                    "updated_widgets": list(updated_widgets),
                }
            )

        if clear_status:
            await self.set_active_cell(None)
        self.cells_with_pending_widget_updates.clear()

        # fixme(rteqs): cleanup signals in some other way. the below does not work because widget signals
        # are restored on `init` but there are no corresponding `rnodes`
        # for x in unused_signals:
        #     del self.widget_signals[x]

    async def update_kernel_snapshot_status(self, key: str, status: str, data: dict[str, int] | None = None) -> None:
        assert key in {"save_kernel_snapshot", "load_kernel_snapshot"}
        self.snapshot_status = status
        msg = {"type": key, "status": status}
        if data is not None:
            msg["data"] = data
        await self.send(msg)

    async def save_kernel_snapshot(self) -> None:
        snapshot_dir.mkdir(parents=True, exist_ok=True)

        s_nodes = {}
        s_signals = {}

        try:
            def add_and_check_listeners(sig: Signal):
                for lid, lis in sig._listeners.items():
                    if lid not in s_nodes:
                        s_nodes[lid] = lis.serialize()
                if sig.id not in s_signals:
                    s_signals[sig.id] = sig.serialize()

            for node in self.cell_rnodes.values():
                s_nodes[node.id] = node.serialize()
                for sig in node.signals.values():
                    add_and_check_listeners(sig)

            for sig in self.widget_signals.values():
                add_and_check_listeners(sig)

            s_globals: dict[str, SerializedSignal | SerializedGlobal] = {}
            for k, val in self.k_globals.items():
                if k in {"__builtins__", "__warningregistry__"}:
                    continue

                if isinstance(val._value, ad.AnnData):
                    s_globals[k] = serialize_anndata(val._value, snapshot_dir)
                elif isinstance(val._value, Signal):
                    s_globals[k] = val._value.serialize()
                elif isinstance(val._value, BaseWidget):
                    if val._value._has_signal:
                        assert val._value._signal.id in s_signals, f"missing {val._value._signal.id}"
                        s_globals[k] = val._value.serialize()
                else:
                    s_val, msg = safe_serialize_obj(val._value)
                    s_globals[k] = {"value": s_val, "error_msg": msg}

            s_anndata = list(self.ann_data_objects.keys())

            s_depens = {
                "s_globals": s_globals,
                "s_nodes": s_nodes,
                "s_signals": s_signals,
                "widget_signals": {k: v.id for k, v in self.widget_signals.items()},
                "nodes_with_widgets": list(self.nodes_with_widgets.keys()),
                "cell_rnodes": {k: v.id for k, v in self.cell_rnodes.items()},
                "s_anndata": s_anndata,
                # todo(kenny): figure out what to do with these
                "ldata_dataframes": {},
                "registry_dataframes": {},
                "url_dataframes": {},
            }

            data = orjson.dumps(s_depens)
        except Exception:
            await self.update_kernel_snapshot_status("save_kernel_snapshot", "error", {"error_msg": traceback.format_exc()})
            return

        total = len(data)

        await self.update_kernel_snapshot_status("save_kernel_snapshot", "start", {"progress_bytes": 0, "total_bytes": total})

        with (snapshot_dir / snapshot_f_name).open("wb") as f:
            saved_since_last = 0
            for start in range(0, total, snapshot_chunk_bytes):
                end = min(start + snapshot_chunk_bytes, total)
                f.write(data[start:end])

                saved_since_last += (end - start)
                if saved_since_last >= snapshot_progress_interval_bytes or end == total:
                    await self.update_kernel_snapshot_status("save_kernel_snapshot", "progress", {"progress_bytes": end, "total_bytes": total})
                    saved_since_last = 0

        await self.update_kernel_snapshot_status("save_kernel_snapshot", "done")

    async def load_kernel_snapshot(self) -> None:

        snapshot_f = snapshot_dir / snapshot_f_name
        if not snapshot_f.exists():
            await self.update_kernel_snapshot_status("load_kernel_snapshot", "done")
            return

        total = snapshot_f.stat().st_size
        await self.update_kernel_snapshot_status("load_kernel_snapshot", "start", {"progress_bytes": 0, "total_bytes": total})

        data = bytearray()
        read_since_last = 0
        with snapshot_f.open("rb") as f:
            while True:
                chunk = f.read(snapshot_chunk_bytes)
                if chunk == b"":
                    break

                data.extend(chunk)
                read_since_last += len(chunk)

                if read_since_last >= snapshot_progress_interval_bytes or len(data) == total:
                    await self.update_kernel_snapshot_status("load_kernel_snapshot", "progress", {"progress_bytes":  len(data), "total_bytes": total})
                    read_since_last = 0

        try:
            s_depens = orjson.loads(data)

            s_nodes: dict[str, SerializedNode] = s_depens["s_nodes"]
            s_signals = s_depens["s_signals"]

            nodes: dict[str, Node] = {}
            signals: dict[str, Signal[object]] = {}

            for nid, s_node in s_nodes.items():
                nodes[nid] = Node.load(s_node)

            for sid, s_sig in s_signals.items():
                signals[sid] = Signal.load(s_sig)

            for nid, s_node in s_nodes.items():
                node = nodes[nid]

                # note(kenny): node children will get constructed by parents when
                # called so no need to build entire reactive tree
                node.signals = {x: signals.get(x) for x in s_node["signals"]}

            for sid, s_signal in s_signals.items():
                signal = signals[sid]
                signal._listeners = {x: nodes.get(x) for x in s_signal["listeners"]}

            self.restored_nodes = nodes
            self.restored_signals = signals

            self.cell_rnodes = filter_items(s_depens["cell_rnodes"], nodes)
            self.widget_signals = filter_items(s_depens["widget_signals"], signals)

            restored_globals = {}
            self.ann_data_objects = {}
            for k, s_v in s_depens["s_globals"].items():
                if "listeners" in s_v:
                    sig = signals.get(s_v["id"])
                    if sig is None:
                        sig = Signal.load(s_v)
                    restored_globals[k] = sig.serialize(short_val=True)
                    self.k_globals._direct_set(k, Signal(sig))
                elif "_is_plots_faas_widget" in s_v:
                    widget = load_widget_helper(s_v, signals)
                    restored_globals[k] = widget.serialize()
                    self.k_globals._direct_set(k, Signal(widget))
                elif "_is_anndata" in s_v:
                    adata_key = s_v["key"]
                    if adata_key in self.ann_data_objects:
                        adata = self.ann_data_objects[adata_key]
                    else:
                        adata = load_anndata(s_v, snapshot_dir)
                        self.ann_data_objects[adata_key] = adata
                    restored_globals[k] = s_v
                    self.k_globals._direct_set(k, Signal(adata))
                else:
                    val, error_msg = safe_unserialize_obj(s_v["value"])
                    if val is unable_to_unserialize_symbol:
                        restored_globals[k] = {
                            "value": small_repr(val),
                            "msg": f"unserializable. stored err: {s_v['error_msg']}, unserial err: {error_msg}",
                        }
                    else:
                        restored_globals[k] = {
                            "value": small_repr(val),
                            "msg": f"stored error: {s_v['error_msg']}",
                        }
                        self.k_globals._direct_set(k, Signal(val))

            self.restored_globals = restored_globals

            nodes_with_widgets = {}
            for x in s_depens["nodes_with_widgets"]:
                node = nodes.get(x)
                if node is not None:
                    nodes_with_widgets[x] = node

            self.nodes_with_widgets = nodes_with_widgets
        except Exception:

            await self.update_kernel_snapshot_status("load_kernel_snapshot", "error", {"error_msg": traceback.format_exc()})
            return

        await self.update_kernel_snapshot_status("load_kernel_snapshot", "done")

    def get_widget_value(self, key: str) -> Signal[object]:
        assert ctx.cur_comp is not None

        key = f"{ctx.cur_comp.name_path()}/{key}"

        if key not in self.widget_signals:
            self.widget_signals[key] = Signal(Nothing.x, name=key)

        return self.widget_signals[key]

    def emit_widget(self, key: str, data: WidgetState) -> None:
        assert ctx.cur_comp is not None
        assert loop is not None

        ctx.cur_comp.widget_states[key] = data
        self.nodes_with_widgets[ctx.cur_comp.id] = ctx.cur_comp

        if (data["type"] == "plot" or data["type"] == "table"):
            loop.create_task(self.send({"type": "cell_value_viewer_init", "key": key, "value_viewer_key": data["value_viewer_key"], "global_key": data["global_key"]}))

        # todo(maximsmol): I don't think this is actually nullable anymore
        cell_id = ctx.cur_comp.cell_id
        if cell_id is not None:
            self.cells_with_pending_widget_updates.add(cell_id)

    def submit_widget_state(self) -> None:
        for s in ctx.updated_signals.values():
            s._apply_updates()

        self.conn.call_fut(
            self.on_tick_finished(ctx.signals_updated_from_code, clear_status=False)
        ).result()

    def on_dispose(self, node: Node) -> None:
        if node.id not in self.nodes_with_widgets:
            return

        if node.cell_id is not None:
            self.cells_with_pending_widget_updates.add(node.cell_id)

        del self.nodes_with_widgets[node.id]

    def _cell_outputs(self) -> CategorizedCellOutputs:
        res = CategorizedCellOutputs()
        for x in self.k_globals.available:
            res.all.append(x)

            val = self.k_globals[x]
            if hasattr(val, "iloc"):
                res.dfs.append(x)

            if isinstance(val, BaseFigure):
                res.figures.append(x)

            if isinstance(val, Figure) or (
                hasattr(val, "figure") and isinstance(val.figure, Figure)
            ):
                res.static_figures.append(x)

        res.all.sort()
        res.dfs.sort()
        res.figures.sort()
        res.static_figures.sort()

        return res

    async def exec(self, *, cell_id: str, code: str, _from_stub: bool = False) -> None:
        filename = f"<cell {cell_id}>"

        try:
            if not _from_stub:
                assert ctx.cur_comp is None
                assert not ctx.in_tx

            self.cell_status[cell_id] = "running"

            comp = self.cell_rnodes.get(cell_id)
            if comp is not None and not _from_stub:
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
                # fixme(maximsmol): a cell should also be considered running if a
                # child reactive node is running
                #
                # the only complication is to tell when it has *finished*
                # running so we can set the status & send results + flush logs
                self.cell_status[cell_id] = "running"

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

                except (KeyboardInterrupt, Exception):
                    self.cell_status[cell_id] = "error"
                    await self.send_cell_result(cell_id)

                finally:
                    sys.stdout.flush()
                    sys.stderr.flush()

            x.__name__ = filename

            await ctx.run(x, _cell_id=cell_id, code=code)

        except (KeyboardInterrupt, Exception):
            self.cell_status[cell_id] = "error"
            await self.send_cell_result(cell_id)

    async def send_cell_result(self, cell_id: str) -> None:
        await self.send_global_updates()

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
                "dataframe_json": {
                    "schema": build_table_schema(res, version=False),
                    "type": "pandas",
                },
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
                        "type": "pandas" if not hasattr(res, "compute") else "dask",
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
                {
                    "type": "output_value",
                    **(id_fields),
                    **(key_fields),
                    "webp": res,
                }
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

    async def upload_ldata(
        self,
        *,
        dst: str,
        key: str,
        viewer_id: str | None,
        cell_id: str | None,
        filename: str | None,
    ) -> None:
        # todo(rteqs): figure out how to get value without going through reactive context
        df = self.k_globals[key]

        pagination_settings = self.lookup_pagination_settings(
            cell_id=cell_id, viewer_id=viewer_id, source_key=key
        )

        if df is None or not isinstance(df, pd.DataFrame):
            raise ValueError(
                "viewer does not exists or variable is not a pandas DataFrame"
            )

        res = filter_and_sort(df=df, pagination_settings=pagination_settings)

        tmp_dir = Path.home() / ".latch" / "plots"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        name = filename if filename is not None else key
        local_path = tmp_dir / f"./{name}.csv"
        res.to_csv(local_path, index=True)

        if not local_path.exists():
            raise RuntimeError("unable to save dataframe to csv")

        # todo(rteqs): stream directly to LData without temp file
        LPath(urljoins(dst, local_path.name)).upload_from(local_path)

    async def accept(self) -> None:
        # print("[kernel] accept")
        msg = await self.conn.recv()
        # print("[kernel] <", msg)

        if msg["type"] == "init":
            self.cell_output_selections = msg["cell_output_selections"]
            self.plot_data_selections = msg["plot_data_selections"]
            self.plot_configs = msg["plot_configs"]
            self.session_snapshot_mode = msg["session_snapshot_mode"]

            if self.session_snapshot_mode:
                await self.load_kernel_snapshot()

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

                    sig = self.widget_signals.get(k)
                    if sig is None:
                        self.widget_signals[k] = Signal(v["value"])
                    else:
                        sig._value = v["value"]

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

        if msg["type"] == "upload_ldata":
            await self.upload_ldata(
                dst=msg["dst"],
                key=msg["key"],
                viewer_id=msg.get("viewer_id"),
                cell_id=msg.get("cell_id"),
                filename=msg.get("filename"),
            )

            await self.send({"type": "upload_ldata"})
            return

        if msg["type"] == "save_kernel_snapshot":
            await self.save_kernel_snapshot()
            return

        if msg["type"] == "h5":
            response = await handle_h5_widget_message(msg, self.send)
            if response is not None:
                await self.send(response)
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
    if sys.platform == "linux":
        from ctypes import CDLL

        libc = CDLL("libc.so.6")
        PR_SET_NAME = 15  # https://github.com/torvalds/linux/blob/2df0c02dab829dd89360d98a8a1abaa026ef5798/include/uapi/linux/prctl.h#L56
        libc.prctl(PR_SET_NAME, b"kernel")

    asyncio.run(main())
