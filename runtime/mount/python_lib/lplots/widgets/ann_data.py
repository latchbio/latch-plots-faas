from dataclasses import dataclass
from typing import Literal, NotRequired

import anndata as ad  # type: ignore  # noqa: PGH003

from .. import _inject
from ..reactive import Signal
from . import _emit, _state


class AnnDataState(_emit.WidgetState[Literal["ann_data"], str]):
    obj_id: str
    readonly: NotRequired[bool]


@dataclass(frozen=True, kw_only=True)
class AnnData:
    _key: str
    _state: AnnDataState
    _signal: Signal[object]

    def _value(self, val: object) -> ad.AnnData | None:
        if not isinstance(val, ad.AnnData):
            return None
        return val

    @property
    def value(self) -> ad.AnnData | None:
        res = self._signal()
        return self._value(res)

    def sample(self) -> ad.AnnData | None:
        res = self._signal.sample()
        return self._value(res)


def w_ann_data(
    *,
    key: str | None = None,
    adata: ad.AnnData,
    readonly: bool = False,
) -> AnnData:
    key = _state.use_state_key(key=key)
    obj_id = str(id(adata))

    _inject.kernel.ann_data_objects[obj_id] = adata

    res = AnnData(
        _key=key,
        _state={
            "type": "ann_data",
            "obj_id": obj_id,
            "readonly": readonly,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
