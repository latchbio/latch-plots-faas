from dataclasses import dataclass
from typing import Literal, NotRequired, TypedDict

from ..reactive import Signal
from ..utils.nothing import Nothing
from . import _emit, _state


class ButtonWidgetState(_emit.WidgetState[Literal["button"], str]):
    label: str
    readonly: bool
    default: bool
    last_clicked: int


@dataclass(frozen=True, kw_only=True)
class ButtonWidget:
    _key: str
    _state: ButtonWidgetState
    _signal: Signal[int]

    @property
    def value(self) -> bool:
        res = self._signal()
        if isinstance(res, int) and res > self._state["last_clicked"]:
            self._state["last_clicked"] = res
            return True
        return False

        if res is Nothing.x or not isinstance(res, int):
            return False

        return res


def w_button(
    *,
    key: str | None = None,
    label: str,
    default: bool = False,
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
            "last_clicked": 0,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
