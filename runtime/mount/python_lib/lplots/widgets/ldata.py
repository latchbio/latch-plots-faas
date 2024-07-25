from dataclasses import dataclass
from typing import Literal, NotRequired

from latch.ldata.path import LPath

from ..reactive import Signal
from . import _emit, _state
from .shared import FormInputAppearance


class LDataPickerState(_emit.WidgetState[Literal["ldata_picker"], str]):
    default: NotRequired[str | None]
    label: str
    readonly: bool
    appearance: NotRequired[FormInputAppearance | None]
    required: bool


@dataclass(frozen=True, kw_only=True)
class LDataPicker:
    _key: str
    _state: LDataPickerState
    _signal: Signal[str]

    @property
    def value(self) -> LPath | None:
        res = self._signal()
        if res is None or not isinstance(res, str) or not res.startswith("latch://"):
            res = self._state.get("default")
            if res is None:
                return None

        return LPath(res)


def w_ldata_picker(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    appearance: FormInputAppearance | None = None,
    default: str | None = None,
    required: bool = False,
) -> LDataPicker:
    key = _state.use_state_key(key=key)

    res = LDataPicker(
        _key=key,
        _state={
            "type": "ldata_picker",
            "readonly": readonly,
            "label": label,
            "default": default,
            "appearance": appearance,
            "required": required,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
