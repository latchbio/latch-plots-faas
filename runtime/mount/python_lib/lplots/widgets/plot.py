from dataclasses import dataclass
from typing import Literal, NotRequired

from matplotlib.axes import Axes
from matplotlib.figure import Figure, SubFigure
from plotly.basedatatypes import BaseFigure
from seaborn import FacetGrid, JointGrid, PairGrid

from .. import _inject
from ..reactive import Signal
from . import _emit, _state, widget

plot_widget_type: Literal["plot"] = "plot"


class PlotState(_emit.WidgetState[plot_widget_type, str]):
    label: NotRequired[str | None]
    plot_title: NotRequired[str | None]
    value_viewer_key: str
    global_key: str


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
    source: Figure | SubFigure | Axes | BaseFigure | FacetGrid | PairGrid | JointGrid,
) -> Plot:
    key = _state.use_state_key(key=key)

    global_key = None
    for k, v in _inject.kernel.k_globals.items():
        if isinstance(v, Signal):
            v = v.sample()

        if id(v) == id(source):
            global_key = k
            break

    plot_title: str | None = None
    if isinstance(source, (Figure, SubFigure)):
        plot_title = source.axes[0].get_title()
    elif isinstance(source, Axes):
        plot_title = source.get_title()
    elif isinstance(source, (FacetGrid, PairGrid, JointGrid)):
        plot_title = source.figure.get_suptitle()
        if plot_title == "":
            plot_title = " | ".join(ax.get_title() for ax in source.figure.axes)
    else:
        plot_title = source.layout.title.text

    if global_key is None:
        raise ValueError("Could not find source in global variables")

    res = Plot(
        _key=key,
        _state={
            "type": plot_widget_type,
            "label": label,
            "plot_title": plot_title,
            "value_viewer_key": f"{global_key}_{key}",
            "global_key": global_key,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
