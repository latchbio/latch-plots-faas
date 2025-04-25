from dataclasses import dataclass
from typing import Literal

from . import _emit, _state
from .interactive_widgets import InteractiveWidget


class ColumnWidgetState(_emit.WidgetState[Literal["column"], None]):
    items: list[str]


@dataclass(frozen=True, kw_only=True)
class Column:
    _key: str
    _state: ColumnWidgetState


def w_column(*, key: str | None = None, items: list[InteractiveWidget]) -> Column:
    key = _state.use_state_key(key=key)

    res = Column(
        _key=key,
        _state={"type": "column", "items": [i._key for i in items]},
    )

    _emit.emit_widget(key, res._state)

    return res
