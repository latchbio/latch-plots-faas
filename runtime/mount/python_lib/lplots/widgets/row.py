from typing import Literal

from runtime.mount.lplots.widgets.checkbox import CheckboxWidget
from runtime.mount.lplots.widgets.dataframe import DataframePicker
from runtime.mount.lplots.widgets.ldata import LDataPicker
from runtime.mount.lplots.widgets.multiselect import MultiSelect
from runtime.mount.lplots.widgets.radio import RadioGroups
from runtime.mount.lplots.widgets.registry import RegistryTablePicker
from runtime.mount.lplots.widgets.select import Select

from . import _emit, _state

AllowedWidgets = (
    CheckboxWidget
    | LDataPicker
    | Select
    | RadioGroups
    | MultiSelect
    | CheckboxWidget
    | DataframePicker
    | RegistryTablePicker
)


class RowWidgetState(_emit.WidgetState[Literal["row"], None]):
    items: list[str]


def w_row(
    *,
    key: str | None = None,
    items: set[AllowedWidgets],
) -> None:
    key = _state.use_state_key(key=key)

    w: RowWidgetState = {
        "type": "row",
        "items": [i._key for i in items],
    }

    _emit.emit_widget(key, w)
