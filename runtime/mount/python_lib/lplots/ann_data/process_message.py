from pathlib import Path
from typing import Any

import anndata as ad  # type: ignore  # noqa: PGH003
from lplots.reactive import ctx

ann_data_ops = ["get_embedding_options", "get_embeddings", "get_obs_options", "get_obs"]

ann_data_object_cache: dict[str, ad.AnnData] = {}


def handle_ann_data_widget_message(
    msg: dict[str, str],
) -> dict[str, Any]:
    if msg["type"] != "ann_data" or "key" not in msg:
        return {
            "type": "ann_data",
            "key": "error",
            "error": "Invalid message",
        }

    widget_key = msg["key"]
    assert ctx.cur_comp is not None

    if widget_key not in ctx.cur_comp.widget_states:
        return {
            "type": "ann_data",
            "key": widget_key,
            "error": "Widget not found",
        }

    widget_state = ctx.cur_comp.widget_states[widget_key]
    if "src" not in widget_state:
        return {
            "type": "ann_data",
            "key": widget_key,
            "error": "Widget state does not contain src",
        }

    path_src = Path(widget_state["src"])
    if not path_src.exists():
        return {
            "type": "ann_data",
            "key": widget_key,
            "error": f"File {path_src} does not exist",
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
            "error": f"Invalid operation: {msg.get('op', '`op` key missing from message')}",
        }

    op = msg["op"]

    if op == "get_embedding_options":
        return {
            "type": "ann_data",
            "key": widget_key,
            "options": list(adata.obsm.keys()),
        }

    if op == "get_embeddings":
        if "embedding" not in msg or msg["embedding"] not in adata.obsm:
            return {
                "type": "ann_data",
                "key": widget_key,
                "error": f"Embedding {msg.get('embedding', '`embedding` key missing from message')} not found",
            }

        embedding = adata.obsm[msg["embedding"]]
        embedding = embedding[:, 0:3].tolist()

        return {
            "type": "ann_data",
            "key": widget_key,
            "embedding": embedding,
        }

    if op == "get_obs_options":
        return {
            "type": "ann_data",
            "key": widget_key,
            "options": list(adata.obs.keys()),
        }

    if op == "get_obs":
        if "obs" not in msg or msg["obs"] not in adata.obs:
            return {
                "type": "ann_data",
                "key": widget_key,
                "error": f"Observation {msg.get('obs', '`obs` key missing from message')} not found",
            }

        obs = adata.obs[msg["obs"]]
        return {
            "type": "ann_data",
            "key": widget_key,
            "obs": obs.to_list(),
        }

    raise ValueError(f"Invalid operation: {op}")
