from dataclasses import dataclass
from typing import Literal, NotRequired, TypedDict

from ..reactive import Signal
from . import _emit, _state, widget


class CheckboxInputAppearance(TypedDict, total=False):
    error_text: str | None
    description: str | None


checkbox_type: Literal["checkbox"] = "checkbox"


class CheckboxWidgetState(_emit.WidgetState[checkbox_type, str]):
    label: str
    readonly: bool
    default: bool
    appearance: NotRequired[CheckboxInputAppearance | None]


@dataclass(frozen=True, kw_only=True)
class CheckboxWidget(widget.BaseWidget):
    _key: str
    _state: CheckboxWidgetState
    _signal: Signal[object | bool]

    def _value(self, val: object) -> bool:
        if not isinstance(val, bool):
            return self._state["default"]

        return val

    @property
    def value(self) -> bool:
        res = self._signal()
        return self._value(res)

    def sample(self) -> bool:
        res = self._signal.sample()
        return self._value(res)


_emit.widget_registry[checkbox_type] = CheckboxWidget


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
            "type": checkbox_type,
            "label": label,
            "default": default,
            "readonly": readonly,
            "appearance": appearance,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
