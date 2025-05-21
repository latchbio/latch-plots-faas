from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from runtime.mount.python_lib.lplots.h5.utils import auto_install
from runtime.mount.python_lib.lplots.h5.utils.persistence import use_anndata_key

from .. import _inject
from ..reactive import Signal
from . import _emit, _state, widget
from .shared import OutputAppearance

ad = auto_install.ad

h5_widget_type: Literal["h5"] = "h5"


@dataclass(frozen=True)
class H5AD:
    type: Literal["h5ad"]
    obj_id: str | None


@dataclass(frozen=True)
class H5Spatial:
    type: Literal["h5spatial"]
    path: Path


class H5State(_emit.WidgetState[h5_widget_type, H5AD | H5Spatial]):
    data: H5AD | H5Spatial
    readonly: bool
    appearance: OutputAppearance | None


@dataclass(frozen=True, kw_only=True)
class H5(widget.BaseWidget):
    _key: str
    _state: H5State
    _signal: Signal[object]

    def _value(self, val: object) -> H5AD | H5Spatial | None:
        if not isinstance(val, str):
            return None

        if val not in _inject.kernel.ann_data_objects:
            return None

        return _inject.kernel.ann_data_objects[val]

    @property
    def value(self) -> H5AD | H5Spatial | None:
        res = self._signal()
        return self._value(res)

    def sample(self) -> H5AD | H5Spatial | None:
        res = self._signal.sample()
        return self._value(res)


def w_h5(
    *,
    key: str | None = None,
    ann_data: ad.AnnData | None = None,
    readonly: bool = False,
    appearance: OutputAppearance | None = None,
) -> H5:
    key = _state.use_state_key(key=key)

    anndata_key: str | None = None
    if ann_data is not None:
        anndata_key = use_anndata_key(ann_data)
        _inject.kernel.ann_data_objects[anndata_key] = ann_data

    res = H5(
        _key=key,
        _state={
            "type": h5_widget_type,
            "data": H5AD(type="h5ad", obj_id=anndata_key),
            "readonly": readonly,
            "appearance": appearance,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
