import types
from abc import ABC
from dataclasses import dataclass
from typing import Literal

from . import _emit, _state
from .checkbox import CheckboxWidget
from .dataframe import DataframePicker
from .datasource import TabularDatasourcePicker
from .ldata import LDataPicker
from .multiselect import MultiSelect
from .plot import Plot
from .radio import RadioGroups
from .registry import RegistryTablePicker
from .select import Select
from .table import Table
from .text import TextInputWidget

grid_widget_type: Literal["grid"] = "grid"

AllowedWidgets = (
    CheckboxWidget
    | LDataPicker
    | Select
    | RadioGroups
    | MultiSelect
    | DataframePicker
    | TextInputWidget
    | RegistryTablePicker
    | TabularDatasourcePicker
    | Plot
    | Table
)


# todo(manske): allow for custom css to define layout
class GridChunk:
    items: list[AllowedWidgets]
    col_span: int
    row_span: int


class GridWidgetState(_emit.WidgetState[Literal["grid"], None]):
    grid_items: list[GridChunk]
    columns: int
    rows: int


@dataclass(frozen=True, kw_only=True)
class Grid(ABC):
    _key: str
    _state: GridWidgetState

    def __enter__(self) -> "Grid":
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: types.TracebackType | None,
    ) -> None:
        pass

    def add(
        self,
        items: AllowedWidgets | set[AllowedWidgets],
        col_span: int = 6,
        row_span: int = 1,
    ) -> None:
        grid_items = GridChunk()
        if isinstance(items, AllowedWidgets):
            grid_items.items = [items]
        else:
            grid_items.items = list(items)

        grid_items.col_span = col_span
        grid_items.row_span = row_span
        return self._state["grid_items"].append(grid_items)


# todo(manske): allow for custom css to define layout
def w_grid(*, key: str | None = None, columns: int = 12, rows: int = 1) -> Grid:
    key = _state.use_state_key(key=key)

    res = Grid(
        _key=key,
        _state={
            "type": grid_widget_type,
            "grid_items": [],
            "columns": columns,
            "rows": rows,
        },
    )

    _emit.emit_widget(key, res._state)
    return res
