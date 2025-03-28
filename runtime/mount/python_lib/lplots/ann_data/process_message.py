from typing import Any

import anndata as ad  # type: ignore  # noqa: PGH003
import numpy as np
from numpy.typing import NDArray
from matplotlib.path import Path

from .. import _inject

ann_data_ops = ["init_data", "get_obsm_options", "get_obsm", "get_obs_options", "get_obs", "get_counts_column", "mutate_obs"]

ann_data_object_cache: dict[str, ad.AnnData] = {}
ann_data_index_cache: dict[str, NDArray[np.int64]] = {}
ann_data_var_index_cache: dict[str, tuple[NDArray[np.str_], NDArray[np.str_] | None]] = {}


MAX_VISUALIZATION_CELLS = 100000
RNG = np.random.default_rng()


def get_obsm(
    obj_id: str,
    adata: ad.AnnData,
    obsm_key: str,
) -> tuple[NDArray[np.float32], NDArray[np.str_]] | tuple[None, None]:
    if obsm_key not in adata.obsm:
        return None, None

    n_cells = adata.n_obs

    # todo(aidan): intelligent downsampling to preserve outliers / information in genera
    if n_cells > MAX_VISUALIZATION_CELLS:
        if obj_id not in ann_data_index_cache:
            idxs = RNG.choice(n_cells, size=MAX_VISUALIZATION_CELLS, replace=False)
            ann_data_index_cache[obj_id] = idxs
        else:
            idxs = ann_data_index_cache[obj_id]
    else:
        idxs = np.arange(n_cells)

    # todo(aidan): allow specifying dimensions to fetch (PCA components is an example where this is useful)
    # todo(aidan): support sparse data?
    obsm = np.asarray(adata.obsm[obsm_key][idxs, :2], dtype=np.float32)  # type: ignore  # noqa: PGH003
    index = np.asarray(adata.obs_names[idxs])

    return obsm, index


def get_obs(
    obj_id: str,
    adata: ad.AnnData,
    obs_key: str,
) -> tuple[NDArray[np.str_], tuple[NDArray[np.str_], NDArray[np.int64]], int]:
    obs = np.asarray(adata.obs[obs_key]).astype(str)
    unique_obs, counts = np.unique(obs, return_counts=True)
    n_cells = adata.n_obs

    # todo(aidan): intelligent downsampling to preserve outliers / information in general
    if n_cells > MAX_VISUALIZATION_CELLS:
        if obj_id not in ann_data_index_cache:
            idxs = RNG.choice(n_cells, size=MAX_VISUALIZATION_CELLS, replace=False)
            ann_data_index_cache[obj_id] = idxs
        else:
            idxs = ann_data_index_cache[obj_id]

        obs = obs[idxs]

        truncated_unique_obs = unique_obs
        if len(unique_obs) > MAX_VISUALIZATION_CELLS:
            sorted_indices = np.argsort(-counts)[:MAX_VISUALIZATION_CELLS]
            truncated_unique_obs = unique_obs[sorted_indices]
            counts = counts[sorted_indices]
    else:
        truncated_unique_obs, counts = np.unique(obs, return_counts=True)

    return obs, (truncated_unique_obs, counts), len(unique_obs)


def get_obs_vector(
    obj_id: str,
    adata: ad.AnnData,
    var_index: str,
) -> NDArray[np.int64]:
    n_cells = adata.n_obs

    if n_cells > MAX_VISUALIZATION_CELLS:
        if obj_id not in ann_data_index_cache:
            idxs = RNG.choice(n_cells, size=MAX_VISUALIZATION_CELLS, replace=False)
            ann_data_index_cache[obj_id] = idxs
        else:
            idxs = ann_data_index_cache[obj_id]
    else:
        idxs = np.arange(n_cells)

    return np.asarray(adata[idxs, var_index].to_df().iloc[:, 0:].values.ravel())


def get_var_index(
    obj_id: str,
    adata: ad.AnnData
) -> tuple[NDArray[np.str_], NDArray[np.str_] | None]:
    if obj_id in ann_data_var_index_cache:
        return ann_data_var_index_cache[obj_id]

    var_index = np.asarray(adata.var_names)

    name_key = None
    for key in adata.var_keys():
        if any(key.lower() in x for x in ["gene_symbols", "gene_names", "gene_ids"]):
            name_key = key
            break

    var_names = None
    if name_key is not None:
        var_names = np.asarray(adata.var[name_key])

    ann_data_var_index_cache[obj_id] = var_index, var_names
    return var_index, var_names


def mutate_obs(
    adata: ad.AnnData,
    obsm_key: str,
    obs_key: str,
    obs_value: str | float,
    lasso_points: list[tuple[int, int]],
) -> None:
    embedding = adata.obsm[obsm_key][:, :2]  # type: ignore  # noqa: PGH003
    polygon = Path(lasso_points)

    mask = polygon.contains_points(embedding)
    adata.obs.loc[mask, obs_key] = obs_value


def handle_ann_data_widget_message(
    msg: dict[str, Any],
) -> dict[str, Any]:
    if msg["type"] != "ann_data" or "key" not in msg or "state" not in msg:
        return {
            "type": "ann_data",
            "key": None,
            "value": {
                "error": "Invalid message -- missing `key` or `state`",
            },
        }

    widget_key = msg["key"]
    widget_state: dict[str, Any] = msg["state"]

    if "obj_id" not in widget_state:
        return {
            "type": "ann_data",
            "key": widget_key,
            "value": {
                "error": "Widget state does not contain obj_id",
            },
        }

    obj_id = widget_state["obj_id"]
    if obj_id not in _inject.kernel.ann_data_objects:
        return {
            "type": "ann_data",
            "key": widget_key,
            "value": {
                "error": f"AnnData object with ID {obj_id} not found in cache",
            },
        }

    adata = _inject.kernel.ann_data_objects[obj_id]

    if "op" not in msg or msg["op"] not in ann_data_ops:
        return {
            "type": "ann_data",
            "key": widget_key,
            "value": {
                "error": f"Invalid operation: {msg.get('op', '`op` key missing from message')}",
            },
        }

    op = msg["op"]

    if op == "init_data":
        init_obsm_key = None
        possible_obsm_keys = adata.obsm_keys()
        for key in possible_obsm_keys:
            if "umap" in key.lower():
                init_obsm_key = key
                break
        if init_obsm_key is None and len(possible_obsm_keys) > 0:
            init_obsm_key = possible_obsm_keys[0]

        init_obs_key = None
        possible_obs_keys = adata.obs_keys()
        for key in possible_obs_keys:
            if "cell" in key.lower() and "type" in key.lower():
                init_obs_key = key
                break
        if init_obs_key is None and len(possible_obs_keys) > 0:
            init_obs_key = possible_obs_keys[0]

        obsm = None
        index = None
        if init_obsm_key is not None:
            obsm, index = get_obsm(obj_id, adata, init_obsm_key)

        obs = None
        unique_obs = None
        nrof_obs = None
        counts = None
        if init_obs_key is not None:
            obs, (unique_obs, counts), nrof_obs = get_obs(obj_id, adata, init_obs_key)

        var_index, var_names = get_var_index(obj_id, adata)

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {
                "data": {
                    # display info
                    "num_obs": adata.n_obs,
                    "num_vars": adata.n_vars,

                    # options
                    "possible_obs_keys": possible_obs_keys,
                    "possible_obsm_keys": possible_obsm_keys,

                    # init state with these
                    "init_obs_key": init_obs_key,
                    "init_obsm_key": init_obsm_key,
                    "init_obsm_values": obsm.tolist() if obsm is not None else None,
                    "init_obsm_index": index.tolist() if index is not None else None,
                    "init_obs_values": obs.tolist() if obs is not None else None,
                    "init_obs_unique_values": unique_obs.tolist() if unique_obs is not None else None,
                    "init_obs_counts": counts.tolist() if counts is not None else None,
                    "init_obs_nrof_values": nrof_obs,

                    # var info
                    "init_var_index": var_index.tolist(),
                    "init_var_names": var_names.tolist() if var_names is not None else None,
                }
            },
        }

    if op == "get_obsm_options":
        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {"data": list(adata.obsm.keys())},
        }

    if op == "get_obsm":
        if "obsm_key" not in msg or msg["obsm_key"] not in adata.obsm:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_key,
                "value": {
                    "error": f"Obsm {msg.get('obsm_key', '`obsm_key` key missing from message')} not found",
                },
            }

        obsm, index = get_obsm(obj_id, adata, msg["obsm_key"])
        if obsm is None or index is None:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_key,
                "value": {"error": "Failed to get obsm"},
            }

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {
                "data": {
                    "fetched_for_key": msg["obsm_key"],
                    "obsm": obsm.tolist(),
                    "index": index.tolist(),
                },
            },
        }

    if op == "get_obs_options":
        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {
                "data": list(adata.obs.keys()),
            },
        }

    if op == "get_obs":
        if "obs_key" not in msg or msg["obs_key"] not in adata.obs:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_key,
                "value": {
                    "error": f"Observation {msg.get('obs_key', '`obs_key` key missing from message')} not found",
                },
            }

        obs, (unique_obs, counts), nrof_obs = get_obs(obj_id, adata, msg["obs_key"])

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {
                "data": {
                    "fetched_for_key": msg["obs_key"],
                    "values": obs.tolist(),
                    "unique_values": unique_obs.tolist(),
                    "counts": counts.tolist(),
                    "nrof_values": nrof_obs
                },
            },
        }

    if op == "get_counts_column":
        if "var_index" not in msg or msg["var_index"] not in adata.var_names:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_key,
                "value": {"error": f"Variable {msg.get('var_index', '`var_index` key missing from message')} not found"},
            }

        gene_column = get_obs_vector(obj_id, adata, msg["var_index"])

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {"data": {
                "fetched_for_key": msg["var_index"],
                "values": gene_column.tolist(),
            }},
        }

    if op == "mutate_obs":
        if "obs_key" not in msg:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_key,
                "value": {"error": "`obs_key` key missing from message"},
            }

        obs_key = msg["obs_key"]
        if obs_key not in adata.obs:
            adata.obs = adata.obs.reindex(columns=[*adata.obs.columns.tolist(), obs_key])

        if "obs_value" in msg and "lasso_points" in msg and "obsm_key" in msg:
            mutate_obs(adata, msg["obsm_key"], obs_key, msg["obs_value"], msg["lasso_points"])

        obs, (unique_obs, counts), nrof_obs = get_obs(obj_id, adata, obs_key)
        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {"data": {
                "obs_key": obs_key,
                "values": obs.tolist(),
                "unique_values": unique_obs.tolist(),
                "counts": counts.tolist(),
                "nrof_values": nrof_obs,
            }},
        }

    raise ValueError(f"Invalid operation: {op}")
