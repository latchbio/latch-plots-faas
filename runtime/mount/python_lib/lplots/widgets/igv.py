from dataclasses import dataclass
from typing import Literal, NotRequired

from latch.ldata.path import LPath

from ..reactive import Signal
from . import _emit, _state, widget

igv_type: Literal["igv"] = "igv"


class IGVState(_emit.WidgetState[igv_type, str]):
    default: NotRequired[str | None]
    readonly: bool
    required: bool


@dataclass(frozen=True, kw_only=True)
class IGV(widget.BaseWidget):
    _key: str
    _state: IGVState
    _signal: Signal[object]

    def _value(self, val: object) -> LPath | None:
        if not isinstance(val, str) or not val.startswith("latch://"):
            val = self._state.get("default")
            if val is None:
                return None

        return LPath(val)

    @property
    def value(self) -> LPath | None:
        res = self._signal()
        return self._value(res)

    def sample(self) -> LPath | None:
        res = self._signal.sample()
        return self._value(res)


_emit.widget_registry[igv_type] = IGV


def w_igv(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    required: bool = False,
) -> IGV:
    key = _state.use_state_key(key=key)

    res = IGV(
        _key=key,
        _state={
            "type": igv_type,
            "label": label,
            "readonly": readonly,
            "required": required,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
