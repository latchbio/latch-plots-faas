from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal, NotRequired

from ..reactive import Signal
from . import _emit, _state, widget
from .shared import FormInputAppearance

radio_group_type: Literal["radio_group"] = "radio_group"


class RadioGroupState(_emit.WidgetState[radio_group_type, str]):
    label: str
    readonly: bool
    options: list[str]
    default: NotRequired[str | None]
    appearance: NotRequired[FormInputAppearance | None]
    direction: Literal["horizontal", "vertical"]
    required: bool


@dataclass(frozen=True, kw_only=True)
class RadioGroups(widget.BaseWidget):
    _key: str
    _state: RadioGroupState
    _signal: Signal[object | str]

    def _value(self, val: object) -> str | None:
        if not isinstance(val, str) or val not in self._state["options"]:
            val = self._state.get("default")
            if val is None:
                return None

        return val

    @property
    def value(self) -> str | None:
        res = self._signal()
        return self._value(res)

    def sample(self) -> str | None:
        res = self._signal.sample()
        return self._value(res)


_emit.widget_registry[radio_group_type] = RadioGroups


def w_radio_group(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    options: Iterable[str],
    default: str | None = None,
    appearance: FormInputAppearance | None = None,
    required: bool = False,
    direction: Literal["horizontal", "vertical"] = "horizontal",
) -> RadioGroups:
    key = _state.use_state_key(key=key)

    res = RadioGroups(
        _key=key,
        _state={
            "type": radio_group_type,
            "label": label,
            "readonly": readonly,
            "options": list(options),
            "default": default,
            "appearance": appearance,
            "required": required,
            "direction": direction,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
