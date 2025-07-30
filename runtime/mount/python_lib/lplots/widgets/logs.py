import types
from dataclasses import dataclass
from typing import Literal, NotRequired

from runtime.mount.lplots.widgets.shared import FormInputAppearance

from . import _emit, _state, widget

logs_display_type: Literal["logs_display"] = "logs_display"


class LogsDisplayState(_emit.WidgetState[Literal["logs_display"], None]):
    appearance: NotRequired[FormInputAppearance | None]
    format: bool


@dataclass(frozen=True, kw_only=True)
class LogsDisplay(widget.BaseWidget):
    _key: str
    _state: LogsDisplayState


def w_logs_display(
    *,
    key: str | None = None,
    format: bool = False,
    appearance: FormInputAppearance | None = None,
) -> LogsDisplay:
    key = _state.use_state_key(key=key)

    res = LogsDisplay(
        _key=key,
        _state={
            "type": "logs_display",
            "format": format,
            "appearance": appearance,
        },
    )

    _emit.emit_widget(key, res._state)
    return res
