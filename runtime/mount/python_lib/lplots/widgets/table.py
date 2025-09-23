from dataclasses import dataclass
from typing import Literal, NotRequired

from pandas import DataFrame

from .. import _inject
from ..reactive import Signal
from . import _emit, _state, widget
from .shared import OutputAppearance

table_widget_type: Literal["table"] = "table"


class TableState(_emit.WidgetState[table_widget_type, str]):
    label: NotRequired[str | None]
    value_viewer_key: str
    global_key: str
    appearance: OutputAppearance | None


@dataclass(frozen=True, kw_only=True)
class Table(widget.BaseWidget):
    _key: str
    _state: TableState
    _signal: Signal[object]

    def _value(self, val: object) -> None:
        # todo(aidan,nathan): resolve to df
        return None

    @property
    def value(self) -> None:
        res = self._signal()
        return self._value(res)

    def sample(self) -> None:
        res = self._signal.sample()
        return self._value(res)


_emit.widget_registry[table_widget_type] = Table


def w_table(
    *,
    key: str | None = None,
    label: str | None = None,
    source: DataFrame | None = None,
    appearance: OutputAppearance | None = None,
) -> Table:
    key = _state.use_state_key(key=key)

    global_key = None
    for k, v in _inject.kernel.k_globals.items():
        if isinstance(v, Signal):
            v = v.sample()

        if id(v) == id(source):
            global_key = k
            break

    res = Table(
        _key=key,
        _state={
            "type": table_widget_type,
            "label": label,
            "value_viewer_key": f"{global_key}_{key}",
            "global_key": str(global_key),
            "appearance": appearance,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
