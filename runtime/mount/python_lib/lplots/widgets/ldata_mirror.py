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
    dir: LPath | str


@dataclass(frozen=True, kw_only=True)
class LDataMirror(widget.BaseWidget):
    _key: str
    _state: LDataMirrorState
    _signal: Signal[object | LPath]

    def _value(
        self,
    ) -> LPath | None:
        dir = self._state.get("dir")

        if isinstance(dir, str) and dir.startswith("latch://"):
            lpath = LPath(dir)
            # todo(manske): create the directory if it doesn't exist?
            if lpath.is_dir():
                return lpath
        elif isinstance(dir, LPath) and dir.is_dir():
            return dir

        return None

    @property
    def value(self) -> LPath | None:
        self._signal()
        return self._value()

    def sample(self) -> LPath | None:
        self._signal.sample()
        return self._value()


_emit.widget_registry[ldata_mirror_type] = LDataMirror


def w_ldata_mirror(
    *,
    key: str | None = None,
    label: str,
    dir: str | LPath,
    readonly: bool = False,
    appearance: FormInputAppearance | None = None,
    default: str | None = None,
    required: bool = False,
) -> LDataMirror:
    key = _state.use_state_key(key=key)

    res = LDataMirror(
        _key=key,
        _state={
            "type": ldata_mirror_type,
            "readonly": readonly,
            "label": label,
            "appearance": appearance,
            "required": required,
            "dir": dir,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
