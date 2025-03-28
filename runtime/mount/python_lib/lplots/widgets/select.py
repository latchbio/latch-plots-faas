from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime
from typing import Literal, NotRequired

from ..reactive import Signal
from . import _emit, _state, widget
from .shared import FormInputAppearance

select_widget_type: Literal["select"] = "select"


class SelectState(_emit.WidgetState[select_widget_type, str]):
    label: str
    readonly: bool
    options: list[str | int | float | bool | datetime]
    default: NotRequired[str | int | float | bool | datetime | None]
    appearance: NotRequired[FormInputAppearance | None]
    required: bool


@dataclass(frozen=True, kw_only=True)
class Select(widget.BaseWidget):
    _key: str
    _state: SelectState
    _signal: Signal[object]

    def _value(self, val: object) -> str | int | float | bool | datetime | None:
        if val not in self._state["options"]:
            val = self._state.get("default")
            if val is None:
                return None

        return val

    @property
    def value(self) -> str | int | float | bool | datetime | None:
        res = self._signal()
        return self._value(res)

    def sample(self) -> str | int | float | bool | datetime | None:
        res = self._signal.sample()
        return self._value(res)


_emit.widget_registry[select_widget_type] = Select


def w_select(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    options: Iterable[str | int | float | bool | datetime],
    default: str | int | float | bool | datetime | None = None,
    appearance: FormInputAppearance | None = None,
    required: bool = False,
) -> Select:
    key = _state.use_state_key(key=key)

    res = Select(
        _key=key,
        _state={
            "type": select_widget_type,
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
