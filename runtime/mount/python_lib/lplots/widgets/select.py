from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Literal, NotRequired

from ..reactive import Signal
from . import _emit, _state
from .shared import FormInputAppearance


class SelectState(_emit.WidgetState[Literal["select"], str]):
    label: str
    readonly: bool
    options: list[str]
    option_values: dict[str, Any]
    default: NotRequired[str | None]
    appearance: NotRequired[FormInputAppearance | None]
    required: bool


@dataclass(frozen=True, kw_only=True)
class Select:
    _key: str
    _state: SelectState
    _signal: Signal[str]

    @property
    def value(self) -> Any | None:
        res = self._signal()

        option_values = self._state["option_values"]
        if res is None or not isinstance(res, str) or res not in option_values:
            res = self._state.get("default")
            if res is None:
                return None

        return option_values[res]


def w_select(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    options: Iterable[Any],
    default: Any | None = None,
    appearance: FormInputAppearance | None = None,
    required: bool = False,
) -> Select:
    key = _state.use_state_key(key=key)

    option_values: dict[str, Any] = {}
    for opt in options:
        assert hasattr(opt, "__str__")
        option_values[str(opt)] = opt

    if default is not None:
        assert hasattr(default, "__str__")

    res = Select(
        _key=key,
        _state={
            "type": "select",
            "label": label,
            "readonly": readonly,
            "options": [str(opt) for opt in options],
            "option_values": option_values,
            "default": str(default) if default is not None else None,
            "appearance": appearance,
            "required": required,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
