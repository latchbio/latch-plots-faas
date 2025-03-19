from pathlib import Path
from typing import Any

import anndata as ad  # type: ignore  # noqa: PGH003

ann_data_ops = ["get_embedding_options", "get_embeddings", "get_obs_options", "get_obs"]

ann_data_object_cache: dict[str, ad.AnnData] = {}


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

        embedding = adata.obsm[msg["embedding"]]
        embedding = embedding[:, 0:3].tolist()

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {
                "data": embedding,
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

        obs = adata.obs[msg["obs"]]
        return {
            "type": "ann_data",
            "op": op,
            "key": widget_key,
            "value": {
                "data": obs.to_list(),
            },
        }

    raise ValueError(f"Invalid operation: {op}")
