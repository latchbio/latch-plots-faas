from pathlib import Path
from typing import Any

import anndata as ad  # type: ignore  # noqa: PGH003
import numpy as np
from numpy.typing import NDArray

ann_data_ops = ["init_data", "get_obsm_options", "get_obsm", "get_obs_options", "get_obs"]

ann_data_object_cache: dict[str, ad.AnnData] = {}
ann_data_index_cache: dict[str, NDArray[np.int64]] = {}
ann_data_var_index_cache: dict[str, tuple[NDArray[np.str_], NDArray[np.str_] | None]] = {}


MAX_VISUALIZATION_CELLS = 100000
RNG = np.random.default_rng()


def get_obsm(
    src: str,
    adata: ad.AnnData,
    obsm_key: str,
) -> tuple[NDArray[np.float32], NDArray[np.str_]] | tuple[None, None]:
    if obsm_key not in adata.obsm:
        return None, None

    n_cells = adata.n_obs

    # todo(aidan): intelligent downsampling to preserve outliers / information in genera
    if n_cells > MAX_VISUALIZATION_CELLS:
        if src not in ann_data_index_cache:
            idxs = RNG.choice(n_cells, size=MAX_VISUALIZATION_CELLS, replace=False)
            ann_data_index_cache[src] = idxs
        else:
            idxs = ann_data_index_cache[src]
    else:
        idxs = np.arange(n_cells)

    # todo(aidan): allow specifying dimensions to fetch (PCA components is an example where this is useful)
    # todo(aidan): support sparse data?
    obsm = np.asarray(adata.obsm[obsm_key][idxs, :2], dtype=np.float32)  # type: ignore  # noqa: PGH003
    index = np.asarray(adata.obs_names[idxs])

    return obsm, index


def get_obs(
    src: str,
    adata: ad.AnnData,
    obs_key: str,
) -> tuple[NDArray[np.str_], tuple[NDArray[np.str_], NDArray[np.int64]], int]:
    obs = np.asarray(adata.obs[obs_key]).astype(str)
    unique_obs, counts = np.unique(obs, return_counts=True)
    n_cells = adata.n_obs

    # todo(aidan): intelligent downsampling to preserve outliers / information in general
    if n_cells > MAX_VISUALIZATION_CELLS:
        if src not in ann_data_index_cache:
            idxs = RNG.choice(n_cells, size=MAX_VISUALIZATION_CELLS, replace=False)
            ann_data_index_cache[src] = idxs
        else:
            idxs = ann_data_index_cache[src]

        obs = obs[idxs]

        truncated_unique_obs = unique_obs
        if len(unique_obs) > MAX_VISUALIZATION_CELLS:
            sorted_indices = np.argsort(-counts)[:MAX_VISUALIZATION_CELLS]
            truncated_unique_obs = unique_obs[sorted_indices]
            counts = counts[sorted_indices]
    else:
        truncated_unique_obs, counts = np.unique(obs, return_counts=True)

    return obs, (truncated_unique_obs, counts), len(unique_obs)


def get_var_index(
    src: str,
    adata: ad.AnnData
) -> tuple[NDArray[np.str_], NDArray[np.str_] | None]:
    if src in ann_data_var_index_cache:
        return ann_data_var_index_cache[src]

    var_index = np.asarray(adata.var_names)

    name_key = None
    for key in adata.var_keys():
        if any(key.lower() in x for x in ["gene_symbols", "gene_names", "gene_ids"]):
            name_key = key
            break

    var_names = None
    if name_key is not None:
        var_names = np.asarray(adata.var[name_key])

    ann_data_var_index_cache[src] = var_index, var_names
    return var_index, var_names


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

    if "src" not in widget_state:
        return {
            "type": "ann_data",
            "key": widget_key,
            "value": {
                "error": "Widget state does not contain src",
            },
        }

    path_src = Path(widget_state["src"])
    if not path_src.exists():
        return {
            "type": "ann_data",
            "key": widget_key,
            "value": {
                "error": f"File {path_src} does not exist",
            },
        }

    if widget_state["src"] in ann_data_object_cache:
        adata = ann_data_object_cache[widget_state["src"]]
    else:
        adata = ad.read_h5ad(path_src, backed="r")
        ann_data_object_cache[widget_state["src"]] = adata

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
            obsm, index = get_obsm(widget_state["src"], adata, init_obsm_key)

        obs = None
        unique_obs = None
        nrof_obs = None
        counts = None
        if init_obs_key is not None:
            obs, (unique_obs, counts), nrof_obs = get_obs(widget_state["src"], adata, init_obs_key)

        var_index, var_names = get_var_index(widget_state["src"], adata)

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
                    "var_index": var_index,
                    "var_names": var_names,
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

        obsm, index = get_obsm(widget_state["src"], adata, msg["obsm_key"])
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

        obs, (unique_obs, counts), nrof_obs = get_obs(widget_state["src"], adata, msg["obs_key"])

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

    if op == "get_var_index":
        var_index, var_names = get_var_index(widget_state["src"], adata)

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {
                "data": {
                    "var_index": var_index,
                    "var_names": var_names,
                },
            },
        }

    raise ValueError(f"Invalid operation: {op}")
