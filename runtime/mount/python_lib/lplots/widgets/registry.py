from dataclasses import dataclass
from typing import Literal, NotRequired, TypedDict

from latch.registry.table import Table

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


registry_table_type: Literal["registry_table"] = "registry_table"


class RegistryTableState(_emit.WidgetState[registry_table_type, str]):
    label: str
    readonly: bool
    default: NotRequired[str | None]
    appearance: NotRequired[FormInputAppearance | None]
    table_id: str


# note(manske): typed dict to allow additional values like selection/active filters
class RegistryTableValue(TypedDict):
    table: Table | None


@dataclass(frozen=True, kw_only=True)
class RegistryTable(widget.BaseWidget):
    _key: str
    _state: RegistryTableState
    _signal: Signal[object | str]

    def _value(self) -> Table | None:
        table_id = self._state.get("table_id")

        return Table(id=table_id)

    @property
    def value(self) -> Table | None:
        self._signal()
        return self._value()

    def sample(self) -> Table | None:
        return self._value()


_emit.widget_registry[registry_table_type] = RegistryTable


def w_registry_table(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    default: str | None = None,
    appearance: FormInputAppearance | None = None,
    table_id: str,
) -> RegistryTable:
    key = _state.use_state_key(key=key)

    res = RegistryTable(
        _key=key,
        _state={
            "type": registry_table_type,
            "label": label,
            "readonly": readonly,
            "default": default,
            "appearance": appearance,
            "table_id": table_id,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
