from dataclasses import dataclass
from typing import Literal, NotRequired, TypedDict

from molviewspec.builder import Root, State
from molviewspec.nodes import validate_state_tree

from . import _emit, _state, widget

molstar_type: Literal["molstar"] = "molstar"


class MolstarOptions(TypedDict, total=False):
    """Options for Molstar widget configuration."""


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
        try:
            state_json = molstarviewspec.model_dump_json()
            validate_state_tree(state_json)
        except Exception as e:
            raise ValueError(f"Invalid molstarviewspec state: {e}") from e
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
