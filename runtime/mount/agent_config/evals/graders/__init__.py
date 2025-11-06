from .base import BinaryGrader, GraderResult
from .cell_typing import CellTypingGrader
from .numeric_tolerance import NumericToleranceGrader
from .label_set_jaccard import LabelSetJaccardGrader
from .distribution_comparison import DistributionComparisonGrader
from .marker_gene_precision_recall import MarkerGenePrecisionRecallGrader
from .marker_gene_separation import MarkerGeneSeparationGrader
from .proportion_consistency import ProportionConsistencyGrader

GRADER_REGISTRY = {
    "cell_typing": CellTypingGrader,
    "numeric_tolerance": NumericToleranceGrader,
    "label_set_jaccard": LabelSetJaccardGrader,
    "distribution_comparison": DistributionComparisonGrader,
    "marker_gene_precision_recall": MarkerGenePrecisionRecallGrader,
    "marker_gene_separation": MarkerGeneSeparationGrader,
    "proportion_consistency": ProportionConsistencyGrader,
}

__all__ = [
    "BinaryGrader",
    "GraderResult",
    "CellTypingGrader",
    "NumericToleranceGrader",
    "LabelSetJaccardGrader",
    "DistributionComparisonGrader",
    "MarkerGenePrecisionRecallGrader",
    "MarkerGeneSeparationGrader",
    "ProportionConsistencyGrader",
    "GRADER_REGISTRY",
]
