from dataclasses import dataclass
from typing import Literal, NotRequired, TypedDict

from molviewspec.builder import Root, State

from ..reactive import Signal
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


class SelectedResidue(TypedDict):
    chainId: str
    authChainId: str
    residueName: str
    residueNumber: int
    authResidueNumber: int
    compId: str
    oneLetterCode: str


class SequenceSegment(TypedDict):
    chainId: str
    residues: list[SelectedResidue]
    sequence: str
    startResidue: int
    endResidue: int


class SequenceSelection(TypedDict):
    segments: list[SequenceSegment]
    fullSequence: str
    structureLabel: NotRequired[str]


class MolstarValue(TypedDict):
    selection: SequenceSelection | None


class MolstarState(_emit.WidgetState[molstar_type, MolstarValue]):
    label: NotRequired[str | None]
    options: MolstarOptions
    molstarviewspec: NotRequired[State | None]


@dataclass(frozen=True, kw_only=True)
class Molstar(widget.BaseWidget):
    _key: str
    _state: MolstarState
    _signal: Signal[object | MolstarValue]

    def _value(self, val: object) -> MolstarValue:
        default = MolstarValue(selection=None)

        if not isinstance(val, dict):
            return default

        selection = val.get("selection")
        if not isinstance(selection, dict):
            return default

        segments = selection.get("segments")
        if not isinstance(segments, list):
            return default

        full_sequence = selection.get("fullSequence")
        if not isinstance(full_sequence, str):
            return default

        structure_label = selection.get("structureLabel")

        res_selection = SequenceSelection(segments=segments, fullSequence=full_sequence)
        if isinstance(structure_label, str):
            res_selection["structureLabel"] = structure_label

        return MolstarValue(selection=res_selection)

    @property
    def value(self) -> MolstarValue:
        res = self._signal()
        return self._value(res)

    def sample(self) -> MolstarValue:
        res = self._signal.sample()
        return self._value(res)


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
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
