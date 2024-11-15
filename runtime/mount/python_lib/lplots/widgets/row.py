from typing import Literal

from . import _emit, _state
from .checkbox import CheckboxWidget
from .dataframe import DataframePicker
from .datasource import TabularDatasourcePicker
from .ldata import LDataPicker
from .multiselect import MultiSelect
from .radio import RadioGroups
from .registry import RegistryTablePicker
from .select import Select
from .text import TextInputWidget

AllowedWidgets = (
    CheckboxWidget
    | LDataPicker
    | Select
    | RadioGroups
    | MultiSelect
    | DataframePicker
    | TextInputWidget
    | CheckboxWidget
    | DataframePicker
    | RegistryTablePicker
    | TabularDatasourcePicker
)


class RowWidgetState(_emit.WidgetState[Literal["row"], None]):
    items: list[str]


def w_row(
    *, key: str | None = None, items: set[AllowedWidgets], stream: bool = False
) -> None:
    key = _state.use_state_key(key=key)

    w: RowWidgetState = {"type": "row", "items": [i._key for i in items]}

    _emit.emit_widget(key, w, stream)
