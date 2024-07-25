from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal, NotRequired

from ..reactive import Signal
from . import _emit, _state
from .shared import FormInputAppearance


class MultiSelectState(_emit.WidgetState[Literal["multi_select"], str]):
    label: str
    readonly: bool
    options: list[str]
    default: NotRequired[list[str] | None]
    appearance: NotRequired[FormInputAppearance | None]
    required: bool


@dataclass(frozen=True, kw_only=True)
class MultiSelect:
    _key: str
    _state: MultiSelectState
    _signal: Signal[list[str]]

    @property
    def value(self) -> list[str] | None:
        res = self._signal()
        if (
            res is None
            or not isinstance(res, list)
            or not all(isinstance(x, str) for x in res)
        ):
            res = self._state.get("default")
            if res is None:
                return None

        return [x for x in res if x in self._state["options"]]


def w_multi_select(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    appearance: FormInputAppearance | None = None,
    options: Iterable[str],
    default: list[str] | None = None,
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
            "default": default,
            "appearance": appearance,
            "required": required,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
