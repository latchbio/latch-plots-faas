from dataclasses import dataclass
from pathlib import Path
from typing import Literal, NotRequired

from ..reactive import Signal
from . import _emit, _state


class AnnDataState(_emit.WidgetState[Literal["ann_data"], str]):
    src: Path
    readonly: NotRequired[bool]


@dataclass(frozen=True, kw_only=True)
class AnnData:
    _key: str
    _state: AnnDataState
    _signal: Signal[object]

    def _value(self, val: object) -> Path | None:
        if not isinstance(val, Path):
            return None

        return val

    @property
    def value(self) -> Path | None:
        res = self._signal()
        return self._value(res)

    def sample(self) -> Path | None:
        res = self._signal.sample()
        return self._value(res)


def w_ann_data(
    *,
    key: str | None = None,
    src: Path,
    readonly: bool = False,
) -> AnnData:
    key = _state.use_state_key(key=key)

    res = AnnData(
        _key=key,
        _state={
            "type": "ann_data",
            "src": src,
            "readonly": readonly,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
