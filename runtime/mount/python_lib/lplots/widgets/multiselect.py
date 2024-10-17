from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Literal, NotRequired

from ..reactive import Signal
from . import _emit, _state
from .shared import FormInputAppearance


class MultiSelectState(_emit.WidgetState[Literal["multi_select"], str]):
    label: str
    readonly: bool
    options: list[str]
    option_values: dict[str, Any]
    default: NotRequired[list[str] | None]
    appearance: NotRequired[FormInputAppearance | None]
    required: bool


@dataclass(frozen=True, kw_only=True)
class MultiSelect:
    _key: str
    _state: MultiSelectState
    _signal: Signal[list[str]]

    @property
    def value(self) -> list[Any] | None:
        res = self._signal()
        option_values = self._state["option_values"]

        if (
            res is None
            or not isinstance(res, list)
            or not all(isinstance(x, str) for x in res)
            or not all(x in option_values for x in res)
        ):
            res = self._state.get("default")
            if res is None:
                return None

        return [option_values[x] for x in res if x in option_values]


def w_multi_select(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    appearance: FormInputAppearance | None = None,
    options: Iterable[Any],
    default: Iterable[Any] | None = None,
    required: bool = False,
) -> MultiSelect:
    key = _state.use_state_key(key=key)

    option_values: dict[str, Any] = {}

    for opt in options:
        assert hasattr(opt, "__str__")
        option_values[str(opt)] = opt

    if default is not None:
        for val in default:
            assert hasattr(val, "__str__")

    res = MultiSelect(
        _key=key,
        _state={
            "type": "multi_select",
            "label": label,
            "readonly": readonly,
            "options": [str(opt) for opt in options],
            "option_values": option_values,
            "default": [str(val) for val in default] if default is not None else None,
            "appearance": appearance,
            "required": required,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
