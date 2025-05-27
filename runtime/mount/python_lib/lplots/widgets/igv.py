from dataclasses import dataclass
from typing import Any, Literal, NotRequired, TypedDict

from . import _emit, _state, widget

igv_type: Literal["igv"] = "igv"


# note(rahul): see https://igv.org/doc/igvjs/#tracks/Tracks
class BaseTrackOptions(TypedDict):
    type: NotRequired[str]
    sourceType: NotRequired[str]
    format: NotRequired[str]
    name: str
    url: NotRequired[str]
    indexURL: NotRequired[str]
    indexed: NotRequired[bool]
    order: NotRequired[int]
    color: NotRequired[str]
    height: NotRequired[int]
    minHeight: NotRequired[int]
    maxHeight: NotRequired[int]
    visibilityWindow: NotRequired[int]
    removable: NotRequired[bool]
    headers: NotRequired[dict[str, str]]
    oathToken: NotRequired[str]


class AnnotationTrackOptions(BaseTrackOptions):
    type: Literal["annotation"]
    displayMode: NotRequired[Literal["SQUISHED", "EXPANDED", "COLLAPSED"]]
    expandedRowHeight: NotRequired[int]
    squishedRowHeight: NotRequired[int]
    nameField: NotRequired[str]
    maxRows: NotRequired[int]
    searchable: NotRequired[bool]
    searchableFields: NotRequired[list[str]]
    filterTypes: NotRequired[list[str]]
    color: NotRequired[str]
    altColor: NotRequired[str]
    colorBy: NotRequired[str]
    colorTable: NotRequired[dict[str, str]]


class AlignmentSortOptions(TypedDict):
    chr: NotRequired[str]
    position: NotRequired[int]
    option: NotRequired[
        Literal["BASE", "STRAND", "INSERT_SIZE", "MATE_CHR", "MQ", "TAG"]
    ]
    tag: NotRequired[str]
    direction: NotRequired[Literal["ASC", "DESC"]]


class AlignmentFilterOptions(TypedDict):
    vendorFailed: NotRequired[bool]
    duplicates: NotRequired[bool]
    secondary: NotRequired[bool]
    supplementary: NotRequired[bool]
    mq: NotRequired[int]
    readgroups: NotRequired[list[str]]


class AlignmentTrackOptions(BaseTrackOptions):
    type: Literal["alignment"]
    showCoverage: NotRequired[bool]
    showAlignment: NotRequired[bool]
    viewAsPairs: NotRequired[bool]
    pairsSupported: NotRequired[bool]
    color: NotRequired[str]
    deletionColor: NotRequired[str]
    skippedColor: NotRequired[str]
    insertionColor: NotRequired[str]
    negStrandColor: NotRequired[str]
    posStrandColor: NotRequired[str]
    pairConnectorColor: NotRequired[str]
    colorBy: NotRequired[str]
    groupBy: NotRequired[str]
    samplingWindowSize: NotRequired[int]
    sampleDepth: NotRequired[int]
    alignmentRowHeight: NotRequired[int]
    readGroup: NotRequired[str]
    sort: NotRequired[AlignmentSortOptions]
    filter: NotRequired[AlignmentFilterOptions]
    showSoftClips: NotRequired[bool]
    showMismatches: NotRequired[bool]
    showAllBases: NotRequired[bool]
    showInsertionText: NotRequired[bool]
    insertionTextColor: NotRequired[str]
    deletionTextColor: NotRequired[str]
    displayMode: NotRequired[Literal["SQUISHED", "EXPANDED", "FULL"]]
    alignmentRowHeight: NotRequired[int]
    squishedRowHeight: NotRequired[int]
    coverageColor: NotRequired[str]
    coverageTrackHeight: NotRequired[int]
    autoscale: NotRequired[bool]
    autoscaleGroup: Any
    min: NotRequired[int]
    max: NotRequired[int]


# todo(rahul): add more track types as needed
TrackOptions = AnnotationTrackOptions | AlignmentTrackOptions | BaseTrackOptions


class ReferenceOptions(TypedDict):
    id: NotRequired[str]
    name: NotRequired[str]
    fastaURL: NotRequired[str]
    indexURL: NotRequired[str]
    compressedIndexURL: NotRequired[str]
    twoBitURL: NotRequired[str]
    cytobandURL: NotRequired[str]
    aliasUrl: NotRequired[str]
    chromSizesURL: NotRequired[str]
    indexed: NotRequired[bool]
    tracks: NotRequired[list[TrackOptions]]
    chromosomeOrder: NotRequired[list[str]]
    headers: NotRequired[dict[str, str]]
    wholeGenome: NotRequired[bool]


class SearchOptions(TypedDict):
    url: NotRequired[str]
    resultsField: NotRequired[str]
    coords: NotRequired[int]
    chromosomeField: NotRequired[str]
    startField: NotRequired[str]
    endField: NotRequired[str]


# todo(rahul): finish properly typing this
class IGVOptions(TypedDict):
    genome: NotRequired[str]
    reference: NotRequired[ReferenceOptions]
    flanking: NotRequired[int]
    genomeList: NotRequired[str]
    loadDefaultGenomes: NotRequired[bool]
    locus: NotRequired[str | list[str]]
    minimumBases: NotRequired[int]
    queryParametersSupported: NotRequired[bool]
    search: NotRequired[SearchOptions]
    showAllChromosomes: NotRequired[bool]
    showChromosomeWidget: NotRequired[bool]
    showNavigation: NotRequired[bool]
    showIdeogram: NotRequired[bool]
    showSVGButton: NotRequired[bool]
    showRuler: NotRequired[bool]
    showCenterGuide: NotRequired[bool]
    showCursorTrackingGuide: NotRequired[bool]
    trackDefaults: NotRequired[dict[str, Any]]
    tracks: NotRequired[list[TrackOptions]]
    roi: NotRequired[list[Any]]
    oauthToken: NotRequired[str]
    apiKey: NotRequired[str]
    clientId: NotRequired[str]
    nucleotideColors: NotRequired[dict[Literal["A", "C", "G", "T", "N"], str]]
    showSampleNames: NotRequired[bool]


class IGVState(_emit.WidgetState[igv_type, str]):
    options: IGVOptions


@dataclass(frozen=True, kw_only=True)
class IGV(widget.BaseWidget):
    _key: str
    _state: IGVState


_emit.widget_registry[igv_type] = IGV


def w_igv(*, key: str | None = None, options: IGVOptions) -> IGV:
    key = _state.use_state_key(key=key)

    res = IGV(_key=key, _state={"type": igv_type, "options": options})
    _emit.emit_widget(key, res._state)

    return res
