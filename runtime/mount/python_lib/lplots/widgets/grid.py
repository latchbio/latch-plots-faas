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

GridChunk = tuple[list[AllowedWidgets], int, int]


class GridWidgetState(_emit.WidgetState[Literal["grid"], None]):
    grid_items: list[GridChunk]
    columns: int
    rows: int


from dataclasses import dataclass


@dataclass(frozen=True, kw_only=True)
class Grid:
    _key: str
    _state: GridWidgetState

    def add(
        self,
        items: AllowedWidgets | set[AllowedWidgets],
        col_span: int = 6,
        row_span: int = 1,
    ) -> None:
        if isinstance(items, AllowedWidgets):
            grid_items: GridChunk = ([items], col_span, row_span)
        else:
            grid_items: GridChunk = (list(items), col_span, row_span)

        return self._state["grid_items"].append(grid_items)


def w_row(*, key: str | None = None, columns: int = 12, rows: int = 1) -> None:
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
