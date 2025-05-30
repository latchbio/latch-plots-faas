from dataclasses import dataclass
from typing import Literal, NotRequired

from ..reactive import Signal
from . import _emit, _state, widget
from .shared import FormInputAppearance

registry_picker_type: Literal["registry_table_picker"] = "registry_table_picker"


class RegistryTablePickerState(_emit.WidgetState[registry_picker_type, str]):
    label: str
    readonly: bool
    default: NotRequired[str | None]
    appearance: NotRequired[FormInputAppearance | None]
    required: bool


@dataclass(frozen=True, kw_only=True)
class RegistryTablePicker(widget.BaseWidget):
    _key: str
    _state: RegistryTablePickerState
    _signal: Signal[object | str]

    def _value(self, val: object) -> str | None:
        if not isinstance(val, str):
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


_emit.widget_registry[registry_picker_type] = RegistryTablePicker


def w_registry_table_picker(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    default: str | None = None,
    appearance: FormInputAppearance | None = None,
    required: bool = False,
) -> RegistryTablePicker:
    key = _state.use_state_key(key=key)

    res = RegistryTablePicker(
        _key=key,
        _state={
            "type": registry_picker_type,
            "label": label,
            "readonly": readonly,
            "default": default,
            "appearance": appearance,
            "required": required,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
