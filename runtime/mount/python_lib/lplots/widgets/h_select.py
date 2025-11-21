from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime
from typing import Literal, NotRequired, TypeAlias, TypedDict

from ..reactive import Signal
from . import _emit, _state, widget
from .shared import FormInputAppearance

h_select_widget_type: Literal["h_select"] = "h_select"


class HSelectOptionCategory(TypedDict):
    title: str
    children: "list[HSelectOption]"


HSelectOption: TypeAlias = str | int | float | bool | datetime | HSelectOptionCategory


class HSelectState(_emit.WidgetState[h_select_widget_type, str]):
    label: str
    readonly: bool
    options: list[HSelectOption]
    default: NotRequired[set[str | int | float | bool | datetime] | None]
    appearance: NotRequired[FormInputAppearance | None]


@dataclass(frozen=True, kw_only=True)
class HSelect(widget.BaseWidget):
    _key: str
    _state: HSelectState
    _signal: Signal[object]

    def _value(self, val: object) -> set[str | int | float | bool | datetime] | None:
        if not isinstance(val, list):
            val = self._state.get("default")
            if val is None:
                return None

        all_options = set()
        q = list(self._state["options"])
        while len(q) > 0:
            cur = q.pop()
            if isinstance(cur, dict):
                q.extend(cur["children"])
            else:
                all_options.add(cur)

        return {x for x in val if x in self._state["options"]}

    @property
    def value(self) -> set[str | int | float | bool | datetime] | None:
        res = self._signal()
        return self._value(res)

    def sample(self) -> set[str | int | float | bool | datetime] | None:
        res = self._signal.sample()
        return self._value(res)


_emit.widget_registry[h_select_widget_type] = HSelect


def w_h_select(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    options: Iterable[HSelectOption],
    default: set[str | int | float | bool | datetime] | None = None,
    appearance: FormInputAppearance | None = None,
) -> HSelect:
    key = _state.use_state_key(key=key)

    res = HSelect(
        _key=key,
        _state={
            "type": h_select_widget_type,
            "label": label,
            "readonly": readonly,
            "options": list(options),
            "default": default,
            "appearance": appearance,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
