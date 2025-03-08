from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal, NotRequired

from ..reactive import Signal
from . import _emit, _state
from .shared import FormInputAppearance


class RadioGroupState(_emit.WidgetState[Literal["radio_group"], str]):
    label: str
    readonly: bool
    options: list[str]
    default: NotRequired[str | None]
    appearance: NotRequired[FormInputAppearance | None]
    direction: Literal["horizontal", "vertical"]
    required: bool


@dataclass(frozen=True, kw_only=True)
class RadioGroups:
    _key: str
    _state: RadioGroupState
    _signal: Signal[str]

    def _value(self, val: str | None) -> str | None:
        if val is None or val not in self._state["options"]:
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
            "type": "radio_group",
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
