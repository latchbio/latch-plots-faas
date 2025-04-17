import uuid
from dataclasses import dataclass
from typing import Literal, NotRequired

from ..reactive import Signal
from . import _emit, _state, widget

table_widget_type: Literal["table"] = "table"


class TableState(_emit.WidgetState[table_widget_type, str]):
    label: str
    connection_key: str


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
    label: str,
) -> Table:
    key = _state.use_state_key(key=key)

    res = Table(
        _key=key,
        _state={
            "type": table_widget_type,
            "label": label,
            "connection_key": uuid.uuid4().hex,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
