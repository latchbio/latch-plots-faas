from math import sqrt
from typing import Any
import numpy as np


def precalc_box(trace: Any):
    if any(
        x in trace
        for x in [
            "q1",
            "median",
            "q3",
            "lowerfence",
            "upperfence",
            "mean",
            "sd",
            "notchspan",
        ]
    ):
        return

    boxpoints = trace.get("boxpoints", "outliers")
    if boxpoints == "all":
        # todo(maximsmol): support subsampling points
        return

    if trace.get("quartilemethod", "linear") != "linear":
        # todo(maximsmol): support other quantile methods
        return

    orientation = trace.get("orientation", "v")
    data_axis = "y" if orientation == "v" else "x"
    index_axis = "x" if orientation == "v" else "y"

    if index_axis in trace or f"{index_axis}0" in trace or f"d{index_axis}" in trace:
        # todo(maximsmol): support multibox traces
        # todo(maximsmol): support generated index axes
        return

    data = np.sort(np.array(trace.get(data_axis, [])))
    # todo(maximsmol): do not go through `to_json` to avoid serializing to Python arrays

    mean = np.mean(data)
    sd = np.std(data, mean=mean) * trace.get("sdmultiple", 1.0)

    q1, median, q3 = np.quantile(data, [0.25, 0.5, 0.75], method="hazen")
    iqr = q3 - q1

    fence_thres = 1.5 * iqr
    lower_fence_thres = q1 - fence_thres
    upper_fence_thres = q3 + fence_thres

    lower_fence_idx = np.searchsorted(data, lower_fence_thres, side="right")
    upper_fence_idx = np.searchsorted(data, upper_fence_thres, side="left")

    lower_fence = data[lower_fence_idx]
    upper_fence = data[upper_fence_idx]

    outliers = np.extract(data < lower_fence | data > upper_fence, data)
    # todo(maximsmol): support subsampling outliers

    notch_span = 1.57 * iqr / sqrt(len(data))

    # >>> Update result
    if boxpoints == "outliers":
        trace[data_axis] = outliers
        trace["boxpoints"] = "all"
    else:
        # todo(maximsmol): unsupported
        ...

    trace["q1"] = q1
    trace["median"] = median
    trace["q3"] = q3
    trace["lowerfence"] = lower_fence
    trace["upperfence"] = upper_fence
    trace["mean"] = mean
    trace["sd"] = sd
    trace["notchspan"] = notch_span
