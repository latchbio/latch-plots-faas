from pathlib import Path
from typing import Any

import anndata as ad  # type: ignore  # noqa: PGH003
import numpy as np

ann_data_ops = ["get_embedding_options", "get_embeddings", "get_obs_options", "get_obs"]

ann_data_object_cache: dict[str, ad.AnnData] = {}
ann_data_index_cache: dict[str, np.ndarray] = {}


MAX_VISUALIZATION_CELLS = 100000
RNG = np.random.default_rng()


async def handle_ann_data_widget_message(  # noqa: RUF029
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

    if op == "get_embedding_options":
        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {
                "data": list(adata.obsm.keys()),
            },
        }

    if op == "get_embeddings":
        if "embedding" not in msg or msg["embedding"] not in adata.obsm:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_key,
                "value": {
                    "error": f"Embedding {msg.get('embedding', '`embedding` key missing from message')} not found",
                },
            }

        embedding = np.asarray(adata.obsm[msg["embedding"]])
        obs_names = np.asarray(adata.obs_names)

        n_cells = adata.n_obs

        # todo(aidan): intelligent downsampling to preserve outliers / information in general
        if n_cells > MAX_VISUALIZATION_CELLS:
            if widget_state["src"] not in ann_data_index_cache:
                idxs = RNG.choice(n_cells, size=MAX_VISUALIZATION_CELLS, replace=False)
                ann_data_index_cache[widget_state["src"]] = idxs
            else:
                idxs = ann_data_index_cache[widget_state["src"]]

            embedding = embedding[idxs, :]
            obs_names = obs_names[idxs]

        embedding_list = embedding.tolist()
        obs_names_list = obs_names.tolist()

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {
                "data": {
                    "embedding": embedding_list,
                    "obs_names": obs_names_list,
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

        obs = np.asarray(adata.obs[msg["obs"]])
        n_cells = adata.n_obs

        # todo(aidan): intelligent downsampling to preserve outliers / information in general
        if n_cells > MAX_VISUALIZATION_CELLS:
            if widget_state["src"] not in ann_data_index_cache:
                idxs = RNG.choice(n_cells, size=MAX_VISUALIZATION_CELLS, replace=False)
                ann_data_index_cache[widget_state["src"]] = idxs
            else:
                idxs = ann_data_index_cache[widget_state["src"]]

            obs = obs[idxs]

            unique_obs = np.unique(obs)
            if len(unique_obs) > MAX_VISUALIZATION_CELLS:
                unique_obs = unique_obs[RNG.choice(len(unique_obs), size=MAX_VISUALIZATION_CELLS, replace=False)]

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {
                "data": {
                    "values": obs,
                    "unique_values": unique_obs,
                },
            },
        }

    raise ValueError(f"Invalid operation: {op}")
