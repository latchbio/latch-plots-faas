from dataclasses import dataclass
from typing import Literal, NotRequired

from latch.ldata.path import LPath

from ..reactive import Signal
from . import _emit, _state, widget

igv_type: Literal["igv"] = "igv"


class IGVState(_emit.WidgetState[igv_type, str]):
    lpath: LPath


@dataclass(frozen=True, kw_only=True)
class IGV(widget.BaseWidget):
    _key: str
    _state: IGVState

_emit.widget_registry[igv_type] = IGV


def w_igv(
    *,
    key: str | None = None,
    lpath: LPath,
) -> IGV:
    key = _state.use_state_key(key=key)

    res = IGV(
        _key=key,
        _state={
            "type": igv_type,
            "lpath": lpath,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
