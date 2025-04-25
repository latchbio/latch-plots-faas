from dataclasses import dataclass
from typing import Literal

from . import _emit, _state
from .interactive_widgets import InteractiveWidget


class RowWidgetState(_emit.WidgetState[Literal["row"], None]):
    items: list[str]


@dataclass(frozen=True, kw_only=True)
class Row:
    _key: str
    _state: RowWidgetState


def w_row(*, key: str | None = None, items: list[InteractiveWidget]) -> Row:
    key = _state.use_state_key(key=key)

    res = Row(
        _key=key,
        _state={"type": "row", "items": [i._key for i in items]},
    )

    _emit.emit_widget(key, res._state)

    return res
