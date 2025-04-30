from dataclasses import dataclass
from typing import Literal

from . import _emit, _state, widget
from .widget import BaseWidget

column_widget_type: Literal["column"] = "column"


class ColumnWidgetState(_emit.WidgetState[Literal["column"], None]):
    items: list[str]


@dataclass(frozen=True, kw_only=True)
class Column(widget.BaseWidget):
    _key: str
    _state: ColumnWidgetState


def w_column(*, key: str | None = None, items: list[BaseWidget]) -> Column:
    key = _state.use_state_key(key=key)

    res = Column(
        _key=key,
        _state={"type": column_widget_type, "items": [i._key for i in items]},
    )

    _emit.emit_widget(key, res._state)

    return res
