from dataclasses import dataclass
from typing import Any, Literal, NotRequired, Required, TypedDict

from . import _emit, _state, widget

igv_type: Literal["igv"] = "igv"


# note(rahul): see https://igv.org/doc/igvjs/#tracks/Tracks
class BaseTrackOptions(TypedDict, total=False):
    type: str
    sourceType: str
    format: str
    name: str
    url: str
    indexURL: str
    indexed: bool
    order: int
    color: str
    height: int
    minHeight: int
    maxHeight: int
    visibilityWindow: int
    removable: bool
    headers: dict[str, str]
    oathToken: str


class AnnotationTrackOptions(BaseTrackOptions, total=False):
    type: Required[Literal["annotation"]]
    displayMode: Literal["SQUISHED", "EXPANDED", "COLLAPSED"]
    expandedRowHeight: int
    squishedRowHeight: int
    nameField: str
    maxRows: int
    searchable: bool
    searchableFields: list[str]
    filterTypes: list[str]
    color: str
    altColor: str
    colorBy: str
    colorTable: dict[str, str]


class AlignmentSortOptions(TypedDict, total=False):
    chr: str
    position: int
    option: Literal["BASE", "STRAND", "INSERT_SIZE", "MATE_CHR", "MQ", "TAG"]
    tag: str
    direction: Literal["ASC", "DESC"]


class AlignmentFilterOptions(TypedDict, total=False):
    vendorFailed: bool
    duplicates: bool
    secondary: bool
    supplementary: bool
    mq: int
    readgroups: list[str]


class AlignmentTrackOptions(BaseTrackOptions, total=False):
    type: Required[Literal["alignment"]]
    showCoverage: bool
    showAlignment: bool
    viewAsPairs: bool
    pairsSupported: bool
    color: str
    deletionColor: str
    skippedColor: str
    insertionColor: str
    negStrandColor: str
    posStrandColor: str
    pairConnectorColor: str
    colorBy: str
    groupBy: str
    samplingWindowSize: int
    sampleDepth: int
    alignmentRowHeight: int
    readGroup: str
    sort: AlignmentSortOptions
    filter: AlignmentFilterOptions
    showSoftClips: bool
    showMismatches: bool
    showAllBases: bool
    showInsertionText: bool
    insertionTextColor: str
    deletionTextColor: str
    displayMode: Literal["SQUISHED", "EXPANDED", "FULL"]
    alignmentRowHeight: int
    squishedRowHeight: int
    coverageColor: str
    coverageTrackHeight: int
    autoscale: bool
    autoscaleGroup: Any
    min: int
    max: int


# todo(rahul): add more track types as needed
TrackOptions = AnnotationTrackOptions | AlignmentTrackOptions | BaseTrackOptions


class ReferenceOptions(TypedDict, total=False):
    id: str
    name: str
    fastaURL: str
    indexURL: str
    compressedIndexURL: str
    twoBitURL: str
    cytobandURL: str
    aliasUrl: str
    chromSizesURL: str
    indexed: bool
    tracks: list[TrackOptions]
    chromosomeOrder: list[str]
    headers: dict[str, str]
    wholeGenome: bool


class SearchOptions(TypedDict, total=False):
    url: str
    resultsField: str
    coords: int
    chromosomeField: str
    startField: str
    endField: str


class IGVOptions(TypedDict, total=False):
    genome: str
    reference: ReferenceOptions
    flanking: int
    genomeList: str
    loadDefaultGenomes: bool
    locus: str | list[str]
    minimumBases: int
    queryParametersSupported: bool
    search: SearchOptions
    showAllChromosomes: bool
    showChromosomeWidget: bool
    showNavigation: bool
    showIdeogram: bool
    showSVGButton: bool
    showRuler: bool
    showCenterGuide: bool
    showCursorTrackingGuide: bool
    trackDefaults: dict[str, Any]
    tracks: list[TrackOptions]
    roi: list[Any]
    oauthToken: str
    apiKey: str
    clientId: str
    nucleotideColors: dict[Literal["A", "C", "G", "T", "N"], str]
    showSampleNames: bool


class IGVState(_emit.WidgetState[igv_type, str]):
    label: NotRequired[str | None]
    options: IGVOptions


@dataclass(frozen=True, kw_only=True)
class IGV(widget.BaseWidget):
    _key: str
    _state: IGVState


_emit.widget_registry[igv_type] = IGV


def w_igv(
    *, key: str | None = None, label: str | None = None, options: IGVOptions
) -> IGV:
    key = _state.use_state_key(key=key)

    res = IGV(
        _key=key,
        _state={
            "type": igv_type,
            "options": options,
            "label": label,
        },
    )
    _emit.emit_widget(key, res._state)

    return res
