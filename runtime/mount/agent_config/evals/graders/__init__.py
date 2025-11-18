from .base import BinaryGrader, GraderResult
from .boolean_check import BooleanCheckGrader
from .cell_typing import CellTypingGrader
from .numeric_tolerance import NumericToleranceGrader
from .label_set_jaccard import LabelSetJaccardGrader
from .distribution_comparison import DistributionComparisonGrader
from .marker_gene_precision_recall import MarkerGenePrecisionRecallGrader
from .marker_gene_separation import MarkerGeneSeparationGrader
from .proportion_consistency import ProportionConsistencyGrader
from .spatial_adjacency_simple import SpatialAdjacencySimpleGrader
from .adjusted_mutual_information import AdjustedMutualInformationGrader

GRADER_REGISTRY = {
    "boolean_check": BooleanCheckGrader,
    "cell_typing": CellTypingGrader,
    "numeric_tolerance": NumericToleranceGrader,
    "label_set_jaccard": LabelSetJaccardGrader,
    "distribution_comparison": DistributionComparisonGrader,
    "marker_gene_precision_recall": MarkerGenePrecisionRecallGrader,
    "marker_gene_separation": MarkerGeneSeparationGrader,
    "proportion_consistency": ProportionConsistencyGrader,
    "spatial_adjacency_simple": SpatialAdjacencySimpleGrader,
    "adjusted_mutual_information": AdjustedMutualInformationGrader,
}

__all__ = [
    "BinaryGrader",
    "GraderResult",
    "BooleanCheckGrader",
    "CellTypingGrader",
    "NumericToleranceGrader",
    "LabelSetJaccardGrader",
    "DistributionComparisonGrader",
    "MarkerGenePrecisionRecallGrader",
    "MarkerGeneSeparationGrader",
    "ProportionConsistencyGrader",
    "SpatialAdjacencySimpleGrader",
    "AdjustedMutualInformationGrader",
    "GRADER_REGISTRY",
]
