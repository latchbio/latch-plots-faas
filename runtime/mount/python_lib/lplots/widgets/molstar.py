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
    chain_id: str
    auth_chain_id: str
    residue_name: str
    residue_number: int
    auth_residue_number: int
    comp_id: str
    one_letter_code: str


class SequenceSegment(TypedDict):
    chain_id: str
    residues: list[SelectedResidue]
    sequence: str
    start_residue: int
    end_residue: int


class SequenceSelection(TypedDict):
    segments: list[SequenceSegment]
    full_sequence: str
    structure_label: NotRequired[str]


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

        segments_raw = selection.get("segments")
        if not isinstance(segments_raw, list):
            return default

        full_sequence = selection.get("fullSequence")
        if not isinstance(full_sequence, str):
            return default

        structure_label = selection.get("structureLabel")

        segments: list[SequenceSegment] = []
        for seg in segments_raw:
            if not isinstance(seg, dict):
                continue

            residues_raw = seg.get("residues")
            if not isinstance(residues_raw, list):
                continue

            residues: list[SelectedResidue] = []
            for res in residues_raw:
                if not isinstance(res, dict):
                    continue

                residues.append(
                    SelectedResidue(
                        chain_id=res.get("chainId", ""),
                        auth_chain_id=res.get("authChainId", ""),
                        residue_name=res.get("residueName", ""),
                        residue_number=res.get("residueNumber", 0),
                        auth_residue_number=res.get("authResidueNumber", 0),
                        comp_id=res.get("compId", ""),
                        one_letter_code=res.get("oneLetterCode", ""),
                    )
                )

            segments.append(
                SequenceSegment(
                    chain_id=seg.get("chainId", ""),
                    residues=residues,
                    sequence=seg.get("sequence", ""),
                    start_residue=seg.get("startResidue", 0),
                    end_residue=seg.get("endResidue", 0),
                )
            )

        res_selection = SequenceSelection(segments=segments, full_sequence=full_sequence)
        if isinstance(structure_label, str):
            res_selection["structure_label"] = structure_label

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
    *,
    key: str | None = None,
    label: str | None = None,
    options: MolstarOptions,
    molstarviewspec_builder: Root | None = None,
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
