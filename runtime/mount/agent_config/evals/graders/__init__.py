from .base import BinaryGrader, GraderResult
from .cell_typing import CellTypingGrader
from .numeric_tolerance import NumericToleranceGrader
from .label_set_jaccard import LabelSetJaccardGrader
from .distribution_comparison import DistributionComparisonGrader

GRADER_REGISTRY = {
    "cell_typing": CellTypingGrader,
    "numeric_tolerance": NumericToleranceGrader,
    "label_set_jaccard": LabelSetJaccardGrader,
    "distribution_comparison": DistributionComparisonGrader,
}

__all__ = [
    "BinaryGrader",
    "GraderResult",
    "CellTypingGrader",
    "NumericToleranceGrader",
    "LabelSetJaccardGrader",
    "DistributionComparisonGrader",
    "GRADER_REGISTRY",
]
