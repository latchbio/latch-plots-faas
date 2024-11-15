from dataclasses import dataclass
from typing import Literal, NotRequired

from ..reactive import Signal
from . import _emit, _state
from .shared import FormInputAppearance


class RegistryTablePickerState(
    _emit.WidgetState[Literal["registry_table_picker"], str]
):
    label: str
    readonly: bool
    default: NotRequired[str | None]
    appearance: NotRequired[FormInputAppearance | None]
    required: bool


@dataclass(frozen=True, kw_only=True)
class RegistryTablePicker:
    _key: str
    _state: RegistryTablePickerState
    _signal: Signal[str]

    @property
    def value(self) -> str | None:
        res = self._signal()
        if res is None or not isinstance(res, str):
            res = self._state.get("default")
            if res is None:
                return None

        return res


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
            "type": "registry_table_picker",
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
