from dataclasses import dataclass, field
from typing import Any, Literal, TypedDict

from ..reactive import Signal
from ..utils.nothing import Nothing
from . import _emit, _state


class ButtonWidgetSignalValue(TypedDict):
    clicked: int
    last_clicked: int


class ButtonWidgetState(_emit.WidgetState[Literal["button"], str]):
    label: str
    readonly: bool
    default: ButtonWidgetSignalValue


def _is_value_widget_signal(data: Any) -> bool:
    if not isinstance(data, dict):
        return False

    required_keys = ButtonWidgetSignalValue.__annotations__.keys()

    if not all(key in data for key in required_keys):
        return False

    if not isinstance(data["clicked"], int):
        return False
    if not isinstance(data["last_clicked"], int):
        return False

    return True


@dataclass(kw_only=True)
class ButtonWidget:
    _key: str
    _state: ButtonWidgetState
    _signal: Signal[ButtonWidgetSignalValue]

    _last_clicked_ref: None | int = field(default=None, repr=False)

    @property
    def value(self) -> bool:
        res = self._signal()

        if res is Nothing.x:
            return False

        if not _is_value_widget_signal(res):
            return False

        if self._last_clicked_ref is None:
            self._last_clicked_ref = res["last_clicked"]

        v = res["clicked"]
        if v > self._last_clicked_ref:
            res["last_clicked"] = v
            return True

        return False


def w_button(
    *,
    key: str | None = None,
    label: str,
    default: ButtonWidgetSignalValue = {
        "clicked": 0,
        "last_clicked": 0,
    },
    readonly: bool = False,
) -> ButtonWidget:
    key = _state.use_state_key(key=key)

    res = ButtonWidget(
        _key=key,
        _state={
            "type": "button",
            "label": label,
            "default": default,
            "readonly": readonly,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
