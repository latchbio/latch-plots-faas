from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime
from typing import Literal, NotRequired

from ..reactive import Signal
from . import _emit, _state
from .shared import FormInputAppearance


class MultiSelectState(_emit.WidgetState[Literal["multi_select"], str]):
    label: str
    readonly: bool
    options: list[str | int | float | bool | datetime]
    default: NotRequired[list[str | int | float | bool | datetime] | None]
    appearance: NotRequired[FormInputAppearance | None]
    required: bool


@dataclass(frozen=True, kw_only=True)
class MultiSelect:
    _key: str
    _state: MultiSelectState
    _signal: Signal[object]

    def _value(self, val: object) -> list[str | int | float | bool | datetime] | None:
        if not isinstance(val, list):
            val = self._state.get("default")
            if val is None:
                return None

        return [x for x in val if x in self._state["options"]]

    @property
    def value(self) -> list[str | int | float | bool | datetime] | None:
        res = self._signal()
        return self._value(res)

    def sample(self) -> list[str | int | float | bool | datetime] | None:
        res = self._signal.sample()
        return self._value(res)


def w_multi_select(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    appearance: FormInputAppearance | None = None,
    options: Iterable[str | int | float | bool | datetime],
    default: Iterable[str | int | float | bool | datetime] | None = None,
    required: bool = False,
) -> MultiSelect:
    key = _state.use_state_key(key=key)

    res = MultiSelect(
        _key=key,
        _state={
            "type": "multi_select",
            "label": label,
            "readonly": readonly,
            "options": list(options),
            "default": list(default) if default is not None else None,
            "appearance": appearance,
            "required": required,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
