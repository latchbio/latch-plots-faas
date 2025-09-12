from math import ceil, pi, sqrt
from typing import Any

import numpy as np

from .precalc_box import precalc_box

gaussian_const = 1 / sqrt(2 * pi)


def gaussian(x):
    return gaussian_const * np.exp(-0.5 * x * x)


# todo(maximsmol): handle non-numeric data
def precalc_violin(trace: Any):
    orientation = trace.get("orientation", "v")
    data_axis = "y" if orientation == "v" else "x"
    index_axis = "x" if orientation == "v" else "y"

    if index_axis in trace:
        # todo(maximsmol): support multibox traces
        print('index_axis', index_axis, 'in trace')
        return

    # note(maximsmol): box precalc will replace this with outliers
    # which we want to happen, but we also need the original data
    # todo(maximsmol): avoid sorting a second time in `precalc_box`
    trace_data = np.sort(np.array(trace.get(data_axis, [])))

    # todo(maximsmol): we don't necessarily need all the data here
    if not precalc_box(trace):
        print('precalc_box returned False')
        return

    trace["density"] = []
    trace["maxKDE"] = []
    trace["count"] = []

    means = trace.get("mean")
    spanmode = trace.get("spanmode", "soft")
    og_trace_span = trace.get("span")

    if spanmode != "manual":
        trace["spanmode"] = "manual"
        trace["span"] = []

    for data_i in range(1):
        data = trace_data

        l = len(data)
        trace["count"].append(l)

        mean = means[data_i] if means is not None else np.mean(data)

        q1 = trace["q1"][data_i]
        q3 = trace["q3"][data_i]
        iqr = q3 - q1

        q0 = data[0]
        q4 = data[-1]

        # >>> bandwidth
        bandwidth = trace.get("bandwidth")
        if isinstance(bandwidth, (list, np.ndarray)):
            bandwidth = bandwidth[data_i]

        if bandwidth is None:
            if q4 - q0 != 0:
                # todo(maximsmol): double check that the ddof is correct
                # pretty sure this is what plotly uses
                ssd = np.std(data, mean=mean, ddof=1)

                silverman = 1.059 * min(ssd, iqr / 1.349) * pow(l, -0.2)
                bandwidth = max(silverman, (q4 - q0) / 100)
            else:
                bandwidth = 0
            trace["bandwidth"] = [bandwidth]

        # >>> span
        span_loose = [q0 - 2 * bandwidth, q4 + 2 * bandwidth]

        trace_span = og_trace_span if og_trace_span is not None else span_loose
        if isinstance(trace_span[0], (list, np.ndarray)):
            trace_span = trace_span[data_i]

        if spanmode != "manual":
            if spanmode == "hard":
                trace_span = [q0, q4]
            elif spanmode == "soft":
                trace_span = [q0 - 2 * bandwidth, q4 + 2 * bandwidth]

            trace["span"].append(trace_span)

        # >>> KDE
        factor = 1 / (l * bandwidth)

        def kernel(x: float):
            return factor * np.sum(gaussian((data - x) / bandwidth))

        span = trace_span[1] - trace_span[0]
        n = ceil(span / (bandwidth / 3))
        step = span / n

        density = []
        maxKDE = 0

        i = 0
        t = trace_span[0]
        while i < n and t < trace_span[1] + step / 2:
            kde = kernel(t)
            density.append({"v": kde, "t": t})
            maxKDE = max(maxKDE, kde)

            i += 1
            t += step

        trace["density"].append(density)
        trace["maxKDE"].append(maxKDE)

    # add dummy val to make trace shown
    # trace[data_axis] = np.asarray([-100])
    trace[index_axis] = np.asarray([-100])