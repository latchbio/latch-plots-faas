from dataclasses import dataclass
from typing import Literal

from latch.ldata.path import LPath

from lplots.h5.utils import auto_install
from lplots.h5.utils.persistence import use_anndata_key

from .. import _inject
from ..reactive import Signal
from . import _emit, _state, widget
from .shared import OutputAppearance

ad = auto_install.ad

h5_widget_type: Literal["h5"] = "h5"


class H5State(_emit.WidgetState[h5_widget_type, str | ad.AnnData | None]):
    obj_id: str | None
    spatial_dir: LPath | None
    readonly: bool
    appearance: OutputAppearance | None


@dataclass(frozen=True, kw_only=True)
class H5(widget.BaseWidget):
    _key: str
    _state: H5State
    _signal: Signal[object | list[tuple[float, float]]]

    def _value(self, val: object) -> list[tuple[float, float]] | None:
        if not isinstance(val, list):
            return None

        for item in val:
            if not (isinstance(item, tuple) and
                   len(item) == 2 and
                   isinstance(item[0], float) and
                   isinstance(item[1], float)):
                return None

        return val

    @property
    def value(self) -> list[tuple[float, float]] | None:
        res = self._signal()
        return self._value(res)

    def sample(self) -> list[tuple[float, float]] | None:
        res = self._signal.sample()
        return self._value(res)


def w_h5(
    *,
    key: str | None = None,
    ann_data: ad.AnnData | None = None,
    spatial_dir: LPath | None = None,
    readonly: bool = False,
    appearance: OutputAppearance | None = None,
) -> H5:
    key = _state.use_state_key(key=key)

    obj_id: str | None = None

    if ann_data is not None:
        anndata_key = use_anndata_key(ann_data)
        _inject.kernel.ann_data_objects[anndata_key] = ann_data
        obj_id = anndata_key

    res = H5(
        _key=key,
        _state={
            "type": h5_widget_type,
            "obj_id": obj_id,
            "spatial_dir": spatial_dir,
            "readonly": readonly,
            "appearance": appearance,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
