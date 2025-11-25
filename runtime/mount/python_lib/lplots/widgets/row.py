from dataclasses import dataclass
from typing import Literal, NotRequired, TypedDict

from . import _emit, _state, widget
from .widget import BaseWidget

row_widget_type: Literal["row"] = "row"


class RowAppearance(TypedDict, total=False):
    compact: bool | None


class RowWidgetState(_emit.WidgetState[row_widget_type, None]):
    items: list[str]
    appearance: NotRequired[RowAppearance | None]


@dataclass(frozen=True, kw_only=True)
class Row(widget.BaseWidget):
    _key: str
    _state: RowWidgetState


def w_row(
    *,
    key: str | None = None,
    items: list[BaseWidget],
    appearance: RowAppearance | None = None,
) -> Row:
    key = _state.use_state_key(key=key)

    res = Row(
        _key=key,
        _state={
            "type": row_widget_type,
            "items": [i._key for i in items],
            "appearance": appearance,
        },
    )

    _emit.emit_widget(key, res._state)

    return res
