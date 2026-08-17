from dataclasses import dataclass, field
from typing import Literal, NotRequired

from matplotlib.axes import Axes
from matplotlib.figure import Figure, SubFigure
from plotly.basedatatypes import BaseFigure
from seaborn import FacetGrid, JointGrid, PairGrid

from .. import _inject
from ..reactive import Signal, ctx
from . import _emit, _state, widget
from .shared import OutputAppearance

plot_widget_type: Literal["plot"] = "plot"


class PlotState(_emit.WidgetState[plot_widget_type, str]):
    label: NotRequired[str | None]
    plot_title: NotRequired[str | None]
    appearance: OutputAppearance | None

    value_viewer_key: str
    global_key: str | None


@dataclass(frozen=True, kw_only=True)
class Plot(widget.BaseWidget):
    _key: str
    _state: PlotState

    _has_signal: bool = field(default=False, init=False)


_emit.widget_registry[plot_widget_type] = Plot


def w_plot(
    *,
    key: str | None = None,
    label: str | None = None,
    source: (
        Figure | SubFigure | Axes | BaseFigure | FacetGrid | PairGrid | JointGrid | None
    ) = None,
    appearance: OutputAppearance | None = None,
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
    if isinstance(source, SubFigure | Figure):
        plot_title = source.axes[0].get_title()
    elif isinstance(source, Axes):
        plot_title = source.get_title()
    elif isinstance(source, BaseFigure):
        plot_title = source.layout.title.text
    elif isinstance(source, FacetGrid | PairGrid | JointGrid):
        plot_title = source.figure.get_suptitle()
        if plot_title == "":
            plot_title = " | ".join(ax.get_title() for ax in source.figure.axes)

    assert ctx.cur_comp is not None
    value_viewer_key = f"{ctx.cur_comp.name_path()}/{key}"
    if global_key is not None:
        # The plot viewer tracks the live global variable by name rather than a
        # snapshot of the figure. If two w_plot widgets point at the same
        # variable, reassigning it (or plotting into it again) overwrites the
        # earlier plot. Guard against this by giving each figure variable a
        # single owning widget and surfacing an error to the agent otherwise.
        widget_identity = value_viewer_key
        owners = _inject.kernel.plot_source_owners
        prev_owner = owners.get(global_key)
        if prev_owner is not None and prev_owner != widget_identity:
            raise ValueError(
                f"w_plot: the figure variable `{global_key}` is already "
                "displayed by another w_plot. Reusing one variable for "
                "multiple plots makes them overwrite each other, because the "
                "plot viewer tracks the live global variable by name. Assign "
                f"this figure to a unique variable (e.g. `{global_key}_1` or "
                f"`{global_key}_umap`) and pass that to w_plot."
            )
        owners[global_key] = widget_identity

        value_viewer_key = f"{global_key}_{key}"

    res = Plot(
        _key=key,
        _state={
            "type": plot_widget_type,
            "label": label,
            "plot_title": plot_title,
            "value_viewer_key": value_viewer_key,
            "global_key": str(global_key) if global_key is not None else None,
            "appearance": appearance,
        },
    )
    emit_state = {**res._state}
    if global_key is None:
        emit_state["source"] = source
    _emit.emit_widget(key, emit_state)

    return res
