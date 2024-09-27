from math import sqrt
from typing import Any
import numpy as np


# todo(maximsmol): handle non-numeric data
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
        return False

    boxpoints = trace.get(
        "boxpoints",
        # note(maximsmol): from violin
        trace.get("points", "outliers"),
    )
    if boxpoints == "all":
        # todo(maximsmol): support subsampling points
        return False

    if trace.get("quartilemethod", "linear") != "linear":
        # todo(maximsmol): support other quantile methods
        return False

    orientation = trace.get("orientation", "v")
    data_axis = "y" if orientation == "v" else "x"
    index_axis = "x" if orientation == "v" else "y"

    if index_axis in trace or (f"{index_axis}0" in trace and f"d{index_axis}" in trace):
        # todo(maximsmol): support multibox traces
        # todo(maximsmol): support generated index axes
        return False

    data = np.sort(np.array(trace.get(data_axis, [])))
    # todo(maximsmol): do not go through `to_json` to avoid serializing to Python arrays

    boxmean = trace.get("boxmean", False)
    sizemode = trace.get("sizemode", "quartiles")
    if boxmean is True or boxmean == "sd" or sizemode == "sd":
        trace.pop("boxmean")

        mean = np.mean(data)
        trace["mean"] = [mean]

        if boxmean == "sd" or sizemode == "sd":
            # todo(maximsmol): should we set different DDOF here? this matches plotly rn i think
            sd = np.std(data, mean=mean) * trace.get("sdmultiple", 1.0)
            trace["sd"] = [sd]

    q1, median, q3 = np.quantile(data, [0.25, 0.5, 0.75], method="hazen")
    trace["q1"] = [q1]
    trace["median"] = [median]
    trace["q3"] = [q3]

    iqr = q3 - q1

    fence_thres = 1.5 * iqr
    lower_fence_thres = q1 - fence_thres
    upper_fence_thres = q3 + fence_thres

    lower_fence_idx = np.searchsorted(data, lower_fence_thres, side="left")
    upper_fence_idx = max(np.searchsorted(data, upper_fence_thres, side="right") - 1, 0)

    lower_fence = data[lower_fence_idx]
    upper_fence = data[upper_fence_idx]

    trace["lowerfence"] = [lower_fence]
    trace["upperfence"] = [upper_fence]

    has_notch = trace.get("notched", "notchwidth" in trace or "notchspan" in trace)
    if has_notch:
        trace.pop("notched")

        notch_span = 1.57 * iqr / sqrt(len(data))
        trace["notchspan"] = [notch_span]

    if boxpoints != "all":
        outliers = np.extract((data < lower_fence) | (data > upper_fence), data)
        # todo(maximsmol): support subsampling outliers

        # todo(maximsmol): fixup `trace.ids`, `trace.text` etc.
        trace[data_axis] = [outliers]

        if "boxpoints" not in trace:
            # note(maximsmol): the default for box plots with precomputed
            # statistics is "all" for some reason
            trace["boxpoints"] = "outliers"

    return True
