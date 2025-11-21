from dataclasses import dataclass
from typing import Literal, NotRequired, TypedDict

from molviewspec.builder import Root, State

from . import _emit, _state, widget

molstar_type: Literal["molstar"] = "molstar"


class LayoutShownByDefault(TypedDict, total=False):
    logs: bool
    sequence_viewer: bool
    right_controls: bool
    left_controls: bool


class LayoutOptions(TypedDict, total=False):
    layout_shown_by_default: LayoutShownByDefault


class MolstarOptions(TypedDict, total=False):
    layout: LayoutOptions


class MolstarState(_emit.WidgetState[molstar_type, str]):
    label: NotRequired[str | None]
    options: MolstarOptions
    molstarviewspec: NotRequired[State | None]


@dataclass(frozen=True, kw_only=True)
class Molstar(widget.BaseWidget):
    _key: str
    _state: MolstarState


_emit.widget_registry[molstar_type] = Molstar


def w_molstar(
    *, key: str | None = None, label: str | None = None, options: MolstarOptions, molstarviewspec_builder: Root | None = None
) -> Molstar:
    key = _state.use_state_key(key=key)

    if molstarviewspec_builder is not None:
        molstarviewspec = molstarviewspec_builder.get_state()
    else:
        molstarviewspec = None

    res = Molstar(
        _key=key,
        _state={
            "type": molstar_type,
            "options": options,
            "label": label,
            "molstarviewspec": molstarviewspec,
        },
    )
    _emit.emit_widget(key, res._state)

    return res
