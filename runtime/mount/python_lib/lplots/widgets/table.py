import uuid
from dataclasses import dataclass
from typing import Literal, NotRequired
from pandas import DataFrame

from ..reactive import Signal
from .. import _inject
from . import _emit, _state, widget

table_widget_type: Literal["table"] = "table"


class TableState(_emit.WidgetState[table_widget_type, str]):
    label: NotRequired[str | None]
    value_viewer_key: str
    found: bool


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
    source: DataFrame,
) -> Table:
    key = _state.use_state_key(key=key)

    found = False
    global_key = None
    for k, v in _inject.kernel.k_globals.items():
        if id(v._value) == id(source):
            global_key = k
            found = True
            break

    res = Table(
        _key=key,
        _state={
            "type": table_widget_type,
            "label": label,
            "value_viewer_key": global_key if global_key is not None else uuid.uuid4().hex,
            "found": found,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
