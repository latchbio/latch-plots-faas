from dataclasses import dataclass
from typing import Literal, NotRequired, TypedDict

from ..reactive import Signal
from ..utils.nothing import Nothing
from . import _emit, _state


class CheckboxInputAppearance(TypedDict, total=False):
    error_text: str | None
    description: str | None


class CheckboxWidgetState(_emit.WidgetState[Literal["checkbox"], str]):
    label: str
    readonly: bool
    default: bool
    appearance: NotRequired[CheckboxInputAppearance | None]


@dataclass(frozen=True, kw_only=True)
class CheckboxWidget:
    _key: str
    _state: CheckboxWidgetState
    _signal: Signal[bool]

    def _value(self, val: bool | Nothing) -> bool:
        if val is Nothing.x:
            return self._state["default"]

        return val

    @property
    def value(self) -> bool:
        res = self._signal()
        return self._value(res)

    def sample(self) -> bool:
        res = self._signal.sample()
        return self._value(res)


def w_checkbox(
    *,
    key: str | None = None,
    label: str,
    default: bool = False,
    readonly: bool = False,
    appearance: CheckboxInputAppearance | None = None,
) -> CheckboxWidget:
    key = _state.use_state_key(key=key)

    res = CheckboxWidget(
        _key=key,
        _state={
            "type": "checkbox",
            "label": label,
            "default": default,
            "readonly": readonly,
            "appearance": appearance,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
