from pathlib import Path
from typing import Any

import anndata as ad  # type: ignore  # noqa: PGH003
import numpy as np
from numpy.typing import NDArray

ann_data_ops = ["init_data", "get_obsm_options", "get_obsm", "get_obs_options", "get_obs"]

ann_data_object_cache: dict[str, ad.AnnData] = {}
ann_data_index_cache: dict[str, NDArray[np.int64]] = {}


MAX_VISUALIZATION_CELLS = 100000
RNG = np.random.default_rng()


def get_obsm(
    src: str,
    adata: ad.AnnData,
    obsm_key: str,
) -> tuple[NDArray[np.float64], NDArray[np.str_]] | tuple[None, None]:
    if obsm_key not in adata.obsm:
        return None, None

    obsm = np.asarray(adata.obsm[obsm_key])
    obs_names = np.asarray(adata.obs_names)

    n_cells = adata.n_obs

    # todo(aidan): intelligent downsampling to preserve information in general
    if n_cells > MAX_VISUALIZATION_CELLS:
        if src not in ann_data_index_cache:
            idxs = RNG.choice(n_cells, size=MAX_VISUALIZATION_CELLS, replace=False)
            ann_data_index_cache[src] = idxs
        else:
            idxs = ann_data_index_cache[src]

        obsm = obsm[idxs, :]
        obs_names = obs_names[idxs]

    return obsm, obs_names


def get_obs(
    src: str,
    adata: ad.AnnData,
    obs_key: str,
) -> tuple[NDArray[np.str_], tuple[NDArray[np.str_], NDArray[np.int64]], int]:
    obs = np.asarray(adata.obs[obs_key])
    n_cells = adata.n_obs

    # todo(aidan): intelligent downsampling to preserve outliers / information in general
    if n_cells > MAX_VISUALIZATION_CELLS:
        if src not in ann_data_index_cache:
            idxs = RNG.choice(n_cells, size=MAX_VISUALIZATION_CELLS, replace=False)
            ann_data_index_cache[src] = idxs
        else:
            idxs = ann_data_index_cache[src]

        obs = obs[idxs]

        unique_obs, counts = np.unique(obs, return_counts=True)
        truncated_unique_obs = unique_obs
        if len(unique_obs) > MAX_VISUALIZATION_CELLS:
            # Keep the most frequent categories up to MAX_VISUALIZATION_CELLS
            sorted_indices = np.argsort(-counts)[:MAX_VISUALIZATION_CELLS]
            truncated_unique_obs = unique_obs[sorted_indices]
            counts = counts[sorted_indices]
    else:
        truncated_unique_obs, counts = np.unique(obs, return_counts=True)

    return obs, (truncated_unique_obs, counts), len(unique_obs)


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
        obs_names = None
        if init_obsm_key is not None:
            obsm, obs_names = get_obsm(widget_state["src"], adata, init_obsm_key)

        obs = None
        unique_obs = None
        nrof_obs = None
        if init_obs_key is not None:
            obs, unique_obs, nrof_obs = get_obs(widget_state["src"], adata, init_obs_key)

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
                    "init_obs_values": obs.tolist() if obs is not None else None,
                    "nrof_unique_obs": nrof_obs,
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
        if "obsm" not in msg or msg["obsm"] not in adata.obsm:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_key,
                "value": {
                    "error": f"Obsm {msg.get('obsm', '`obsm` key missing from message')} not found",
                },
            }

        obsm, obs_names = get_obsm(widget_state["src"], adata, msg["obsm"])
        if obsm is None or obs_names is None:
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
                    "obsm": obsm.tolist(),
                    "obs_names": obs_names.tolist(),
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
        if "obs" not in msg or msg["obs"] not in adata.obs:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_key,
                "value": {
                    "error": f"Observation {msg.get('obs', '`obs` key missing from message')} not found",
                },
            }

        obs, (unique_obs, counts), nrof_obs = get_obs(widget_state["src"], adata, msg["obs"])

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {
                "data": {
                    "values": obs.tolist(),
                    "unique_values": unique_obs.tolist(),
                    "counts": counts.tolist(),
                    "nrof_values": nrof_obs
                },
            },
        }

    raise ValueError(f"Invalid operation: {op}")
