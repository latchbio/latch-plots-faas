from dataclasses import dataclass
from typing import Literal, NotRequired

from latch.ldata.path import LPath

from ..reactive import Signal
from . import _emit, _state, widget
from .shared import FormInputAppearance

ldata_mirror_type: Literal["ldata_mirror"] = "ldata_mirror"


class LDataMirrorState(_emit.WidgetState[ldata_mirror_type, str]):
    label: str
    readonly: bool
    appearance: NotRequired[FormInputAppearance | None]
    required: bool
    dir_node_id: str


@dataclass(frozen=True, kw_only=True)
class LDataMirror(widget.BaseWidget):
    _key: str
    _state: LDataMirrorState
    _dir: LPath
    _signal: Signal[object | LPath]

    @property
    def value(self) -> LPath | None:
        self._signal()
        return self._dir

    def sample(self) -> LPath | None:
        return self._dir


_emit.widget_registry[ldata_mirror_type] = LDataMirror


def w_ldata_mirror(
    *,
    key: str | None = None,
    label: str,
    dir: str | LPath,
    readonly: bool = False,
    appearance: FormInputAppearance | None = None,
    required: bool = False,
) -> LDataMirror:
    key = _state.use_state_key(key=key)

    dir_lpath: LPath | None = None
    if isinstance(dir, str):
        dir_lpath = LPath(dir)
    elif isinstance(dir, LPath):
        dir_lpath = dir

    if dir_lpath is None:
        raise ValueError("Invalid argument `dir`: not an LPath object ")

    if not dir_lpath.is_dir:
        raise ValueError("Invalid argument `dir`: LPath object is not a directory")

    dir_node_id = dir_lpath.node_id()
    if dir_node_id is None:
        raise ValueError("Unable to get node id for directory")

    res = LDataMirror(
        _key=key,
        _dir=dir_lpath,
        _state={
            "type": ldata_mirror_type,
            "readonly": readonly,
            "label": label,
            "appearance": appearance,
            "required": required,
            "dir_node_id": dir_node_id,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
