from dataclasses import dataclass
from typing import Literal, TypedDict

from latch.ldata.path import LPath

from lplots.h5.utils import auto_install
from lplots.h5.utils.persistence import use_anndata_key

from .. import _inject
from ..reactive import Signal
from . import _emit, _state, widget
from .shared import OutputAppearance

ad = auto_install.ad

h5_widget_type: Literal["h5"] = "h5"


class ColorByObs(TypedDict):
    type: Literal["obs"]
    key: str


class ColorByVar(TypedDict):
    type: Literal["var"]
    keys: list[str]


class ViewerPreset(TypedDict):
    default_genes: list[str] | None
    default_color_by: ColorByObs | ColorByVar | None
    default_cell_marker_size: int | None


class H5State(_emit.WidgetState[h5_widget_type, str | ad.AnnData | None]):
    obj_id: str | None
    spatial_dir: LPath | None
    readonly: bool
    appearance: OutputAppearance | None
    viewer_presets: ViewerPreset | None


# note(aidan): typed dict to allow additional values
class H5Value(TypedDict):
    lasso_points: list[list[tuple[float, float]]] | None
    lasso_points_obsm: str | None


@dataclass(frozen=True, kw_only=True)
class H5(widget.BaseWidget):
    _key: str
    _state: H5State
    _signal: Signal[object | H5Value]

    def _value(self, val: object) -> H5Value:
        if not isinstance(val, dict):
            return H5Value(lasso_points=None, lasso_points_obsm=None)

        if "lasso_points" not in val:
            return H5Value(lasso_points=None, lasso_points_obsm=None)

        lasso_points = val["lasso_points"]

        if not isinstance(lasso_points, list):
            return H5Value(lasso_points=None, lasso_points_obsm=None)

        for item in lasso_points:
            if not isinstance(item, list):
                return H5Value(lasso_points=None, lasso_points_obsm=None)

            for point in item:
                if not ((isinstance(point, (tuple, list))) and
                        len(point) == 2 and
                        isinstance(point[0], (float, int)) and
                        isinstance(point[1], (float, int))):
                    return H5Value(lasso_points=None, lasso_points_obsm=None)

        return H5Value(lasso_points=lasso_points, lasso_points_obsm=val.get("lasso_points_obsm", None))

    @property
    def value(self) -> H5Value:
        res = self._signal()
        return self._value(res)

    def sample(self) -> H5Value:
        res = self._signal.sample()
        return self._value(res)


def w_h5(
    *,
    key: str | None = None,
    ann_data: ad.AnnData | None = None,
    spatial_dir: LPath | None = None,
    ann_tiles: LPath | None = None,
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
            "ann_tiles": ann_tiles,
            "readonly": readonly,
            "appearance": appearance,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
