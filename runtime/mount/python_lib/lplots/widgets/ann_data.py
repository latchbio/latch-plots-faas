from dataclasses import dataclass
from typing import Literal

from lplots.h5.utils import auto_install
from lplots.h5.utils.persistence import use_anndata_key

from .. import _inject
from ..reactive import Signal
from . import _emit, _state, widget
from .shared import OutputAppearance

ad = auto_install.ad

ann_data_widget_type: Literal["ann_data"] = "ann_data"


class AnnDataState(_emit.WidgetState[ann_data_widget_type, str]):
    obj_id: str | None
    readonly: bool
    appearance: OutputAppearance | None


@dataclass(frozen=True, kw_only=True)
class AnnData(widget.BaseWidget):
    _key: str
    _state: AnnDataState
    _signal: Signal[object]

    def _value(self, val: object) -> ad.AnnData | None:
        if not isinstance(val, str):
            return None

        if val not in _inject.kernel.ann_data_objects:
            return None

        return _inject.kernel.ann_data_objects[val]

    @property
    def value(self) -> ad.AnnData | None:
        res = self._signal()
        return self._value(res)

    def sample(self) -> ad.AnnData | None:
        res = self._signal.sample()
        return self._value(res)


_emit.widget_registry[ann_data_widget_type] = AnnData


def w_ann_data(
    *,
    key: str | None = None,
    ann_data: ad.AnnData | None = None,
    readonly: bool = False,
    appearance: OutputAppearance | None = None,
) -> AnnData:
    key = _state.use_state_key(key=key)

    anndata_key: str | None = None
    if ann_data is not None:
        anndata_key = use_anndata_key(ann_data)
        _inject.kernel.ann_data_objects[anndata_key] = ann_data

    res = AnnData(
        _key=key,
        _state={
            "type": ann_data_widget_type,
            "obj_id": anndata_key,
            "readonly": readonly,
            "appearance": appearance,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
