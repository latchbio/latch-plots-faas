
import uuid
from dataclasses import dataclass
from typing import Literal, NotRequired

from matplotlib.figure import Figure
from plotly.basedatatypes import BaseFigure

from .. import _inject
from ..reactive import Signal
from . import _emit, _state, widget

plot_widget_type: Literal["plot"] = "plot"


class PlotState(_emit.WidgetState[plot_widget_type, str]):
    label: NotRequired[str | None]
    value_viewer_key: str
    found: bool


@dataclass(frozen=True, kw_only=True)
class Plot(widget.BaseWidget):
    _key: str
    _state: PlotState
    _signal: Signal[object]

    def _value(self, val: object) -> None:
        # todo(aidan,nathan): resolve to plot data
        return None

    @property
    def value(self) -> None:
        res = self._signal()
        return self._value(res)

    def sample(self) -> None:
        res = self._signal.sample()
        return self._value(res)


_emit.widget_registry[plot_widget_type] = Plot


def w_plot(
    *,
    key: str | None = None,
    label: str | None = None,
    source: Figure | BaseFigure,
) -> Plot:
    key = _state.use_state_key(key=key)

    found = False
    global_key = None
    for k, v in _inject.kernel.k_globals.items():
        if v == source:
            global_key = k
            found = True
            break

    res = Plot(
        _key=key,
        _state={
            "type": plot_widget_type,
            "label": label,
            "value_viewer_key": global_key if global_key is not None else uuid.uuid4().hex,
            "found": found,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
