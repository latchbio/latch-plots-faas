from dataclasses import dataclass
from typing import Literal, NotRequired

from latch.ldata.path import LPath

from ..reactive import Signal
from . import _emit, _state, widget
from .shared import FormInputAppearance

ldata_picker_type: Literal["ldata_picker"] = "ldata_picker"


class LDataPickerState(_emit.WidgetState[ldata_picker_type, str]):
    default: NotRequired[str | None]
    label: str
    readonly: bool
    appearance: NotRequired[FormInputAppearance | None]
    required: bool
    file_type: NotRequired[Literal["file", "dir", "any"]]


@dataclass(frozen=True, kw_only=True)
class LDataPicker(widget.BaseWidget):
    _key: str
    _state: LDataPickerState
    _signal: Signal[object | LPath]

    def _value(self, val: object) -> LPath | None:
        if not isinstance(val, str) or not val.startswith("latch://"):
            val = self._state.get("default")
            if val is None:
                return None

        res = LPath(val)

        if self._state.get("file_type") == "dir":
            if not res.is_dir():
                return None

        elif self._state.get("file_type") == "file":
            if res.is_dir():
                return None

        return res

    @property
    def value(self) -> LPath | None:
        res = self._signal()
        return self._value(res)

    def sample(self) -> LPath | None:
        res = self._signal.sample()
        return self._value(res)


_emit.widget_registry[ldata_picker_type] = LDataPicker


def w_ldata_picker(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    appearance: FormInputAppearance | None = None,
    default: str | None = None,
    required: bool = False,
    file_type: Literal["file", "dir", "any"] = "any",
) -> LDataPicker:
    key = _state.use_state_key(key=key)

    res = LDataPicker(
        _key=key,
        _state={
            "type": ldata_picker_type,
            "readonly": readonly,
            "label": label,
            "default": default,
            "appearance": appearance,
            "required": required,
            "file_type": file_type,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
