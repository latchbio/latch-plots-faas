from collections.abc import Awaitable, Callable
from typing import Any

from lplots.h5.h5ad.ops import (
    fetch_and_process_image,
    get_obs,
    get_obs_vector,
    get_obsm,
    get_var_index,
    mutate_obs_by_lasso,
    mutate_obs_by_value,
    pil_image_cache,
)
from lplots.h5.utils import auto_install
from lplots.h5.utils.align import align_image

ad = auto_install.ad

alignment_is_running = False


async def process_h5ad_request(
    msg: dict[str, Any],
    widget_session_key: str,
    adata: ad.AnnData,
    obj_id: str,
    send: Callable[[object], Awaitable[None]],
) -> dict[str, Any]:

    global alignment_is_running

    if "op" not in msg or msg["op"] not in {
        "init_data",
        "get_obsm_options",
        "get_obsm",
        "get_obs_options",
        "get_obs",
        "get_counts_column",
        "mutate_obs",
        "drop_obs",
        "rename_obs",
        "fetch_and_process_image",
        "align_image",
        "store_views",
    }:
        return {
            "type": "h5",
            "key": widget_session_key,
            "value": {
                "error": (
                    f"Invalid operation: {msg.get('op', '`op` key missing from message')}"
                ),
            },
        }

    op = msg["op"]

    max_visualization_cells = msg.get("max_visualization_cells", 100000)

    if op == "init_data":
        init_obsm_key = msg.get("obsm_key")
        possible_obsm_keys = adata.obsm_keys()
        if init_obsm_key is None:
            for key in possible_obsm_keys:
                if "umap" in key.lower() or "spatial" in key.lower():
                    init_obsm_key = key
                    break
        if init_obsm_key is None and len(possible_obsm_keys) > 0:
            init_obsm_key = possible_obsm_keys[0]

        possible_obs_keys = adata.obs_keys()

        init_obs_key = msg.get("obs_key")
        init_var_key = msg.get("var_key")
        if init_obs_key is None and init_var_key is None:
            for key in possible_obs_keys:
                if "cell" in key.lower() and "type" in key.lower():
                    init_obs_key = key
                    break

        if init_obs_key is None and init_var_key is None and len(possible_obs_keys) > 0:
            init_obs_key = possible_obs_keys[0]

        obsm = None
        index = None
        recomputed_index = False
        filters = None
        if init_obsm_key is not None:
            filters = msg.get("filters")
            obsm, index, recomputed_index = get_obsm(
                obj_id, adata, init_obsm_key, filters, max_visualization_cells
            )

        obs = None
        unique_obs = None
        nrof_obs = None
        counts = None
        if init_obs_key is not None and init_obs_key in adata.obs:
            obs, (unique_obs, counts), nrof_obs = get_obs(
                obj_id, adata, init_obs_key, max_visualization_cells
            )

        gene_column = None
        if (
            init_var_key is not None
            and init_obs_key is None
            and init_var_key in adata.var_names
        ):
            gene_column = get_obs_vector(obj_id, adata, init_var_key)

        var_index, var_names = get_var_index(obj_id, adata)

        global alignment_is_running

        return {
            "type": "h5",
            "op": op,
            "data_type": "h5ad",
            "key": widget_session_key,
            "value": {
                "data": {
                    # display info
                    "num_obs": adata.n_obs,
                    "num_vars": adata.n_vars,
                    # options
                    "possible_obs_keys": possible_obs_keys,
                    "possible_obs_keys_types": [
                        str(adata.obs[key].dtype) for key in possible_obs_keys
                    ],
                    "possible_obsm_keys": possible_obsm_keys,
                    # init state with these
                    "init_obs_key": init_obs_key,
                    "init_obsm_key": init_obsm_key,
                    "init_recomputed_index": recomputed_index,
                    "init_obsm_values": obsm.tolist() if obsm is not None else None,
                    "init_obsm_index": index.tolist() if index is not None else None,
                    "init_obsm_filters": filters,
                    "init_obs_values": obs.tolist() if obs is not None else None,
                    "init_obs_unique_values": (
                        unique_obs.tolist() if unique_obs is not None else None
                    ),
                    "init_obs_counts": counts.tolist() if counts is not None else None,
                    "init_obs_nrof_values": nrof_obs,
                    # var info
                    "init_var_index": var_index.tolist(),
                    "init_var_names": (
                        var_names.tolist() if var_names is not None else None
                    ),
                    # var color by info
                    "init_var_values": (
                        gene_column.tolist() if gene_column is not None else None
                    ),
                    "init_var_key": init_var_key if init_var_key is not None else None,
                    # alignment info
                    "alignment_is_running": alignment_is_running,
                    # views info
                    "init_views": adata.uns.get("latch_views", []),
                }
            },
        }

    if op == "get_obsm_options":
        return {
            "type": "h5",
            "op": op,
            "data_type": "h5ad",
            "key": widget_session_key,
            "value": {"data": list(adata.obsm.keys())},
        }

    if op == "get_obsm":
        if "obsm_key" not in msg or msg["obsm_key"] not in adata.obsm:
            return {
                "type": "h5",
                "op": op,
                "data_type": "h5ad",
                "key": widget_session_key,
                "value": {
                    "error": (
                        "Obsm not found"
                        if "obsm_key" in msg
                        else "`obsm_key` key missing from message"
                    )
                },
            }

        filters = msg.get("filters")
        obsm, index, recomputed_index = get_obsm(
            obj_id, adata, msg["obsm_key"], filters, max_visualization_cells
        )
        if obsm is None or index is None:
            return {
                "type": "h5",
                "op": op,
                "data_type": "h5ad",
                "key": widget_session_key,
                "value": {"error": "Failed to get obsm"},
            }

        obs = None
        fetched_for_obs_key = None
        unique_obs = None
        counts = None
        nrof_obs = None
        gene_columns: list[list[str]] | None = None
        fetched_for_var_keys: list[str] | None = None
        if "colored_by_type" in msg and "colored_by_key" in msg:
            if msg["colored_by_type"] == "obs" and msg["colored_by_key"] in adata.obs:
                obs, (unique_obs, counts), nrof_obs = get_obs(
                    obj_id, adata, msg["colored_by_key"], max_visualization_cells
                )
                fetched_for_obs_key = msg["colored_by_key"]
            elif (
                msg["colored_by_type"] == "var"
                and msg["colored_by_key"] in adata.var_names
            ):
                gene_columns = [get_obs_vector(obj_id, adata, msg["colored_by_key"])]
                fetched_for_var_keys = [msg["colored_by_key"]]
            elif (
                msg["colored_by_type"] == "multi_gene"
                and isinstance(msg["colored_by_key"], (list, tuple))
            ):
                fetched_for_var_keys = [
                        g for g in msg["colored_by_key"] if g in
                        adata.var_names
                ]
                gene_columns = [
                    get_obs_vector(obj_id, adata, g) for g in fetched_for_var_keys
                ]

        return {
            "type": "h5",
            "op": op,
            "data_type": "h5ad",
            "key": widget_session_key,
            "value": {
                "data": {
                    "fetched_for_key": msg["obsm_key"],
                    "obsm": obsm.tolist(),
                    "index": index.tolist(),
                    "recomputed_index": recomputed_index,
                    "fetched_for_obs_key": fetched_for_obs_key,
                    "fetched_for_var_keys": fetched_for_var_keys,
                    "fetched_for_filters": filters,
                    "values": obs.tolist() if obs is not None else None,
                    "unique_values": (
                        unique_obs.tolist() if unique_obs is not None else None
                    ),
                    "counts": counts.tolist() if counts is not None else None,
                    "nrof_values": nrof_obs,
                    "var_values": [
                        col.tolist() for col in gene_columns
                    ] if gene_columns is not None else None,
                },
            },
        }

    if op == "get_obs_options":
        return {
            "type": "h5",
            "op": op,
            "data_type": "h5ad",
            "key": widget_session_key,
            "value": {
                "data": list(adata.obs.keys()),
            },
        }

    if op == "get_obs":
        if "obs_key" not in msg or msg["obs_key"] not in adata.obs:
            return {
                "type": "h5",
                "op": op,
                "data_type": "h5ad",
                "key": widget_session_key,
                "value": {
                    "error": (
                        "Observation not found"
                        if "obs_key" in msg
                        else "`obs_key` key missing from message"
                    )
                },
            }

        obs, (unique_obs, counts), nrof_obs = get_obs(
            obj_id, adata, msg["obs_key"], max_visualization_cells
        )

        return {
            "type": "h5",
            "op": op,
            "data_type": "h5ad",
            "key": widget_session_key,
            "value": {
                "data": {
                    "fetched_for_key": msg["obs_key"],
                    "values": obs.tolist(),
                    "unique_values": unique_obs.tolist(),
                    "counts": counts.tolist(),
                    "nrof_values": nrof_obs,
                },
            },
        }

    if op == "get_counts_column":
        if "var_index" not in msg or msg["var_index"] not in adata.var_names:
            return {
                "type": "h5",
                "op": op,
                "key": widget_session_key,
                "value": {
                    "error": (
                        "Variable not found"
                        if "var_index" in msg
                        else "`var_index` key missing from message"
                    )
                },
            }

        gene_column = get_obs_vector(obj_id, adata, msg["var_index"])

        return {
            "type": "h5",
            "op": op,
            "key": widget_session_key,
            "value": {
                "data": {
                    "fetched_for_key": msg["var_index"],
                    "values": gene_column.tolist(),
                }
            },
        }

    if op == "mutate_obs":
        if "obs_key" not in msg:
            return {
                "type": "h5",
                "op": op,
                "data_type": "h5ad",
                "key": widget_session_key,
                "value": {"error": "`obs_key` key missing from message"},
            }

        obs_key = str(msg["obs_key"])
        created_for_key = None
        if obs_key not in adata.obs:
            obs_dtype = msg.get("obs_dtype")
            if obs_dtype is not None and obs_dtype not in {
                "category",
                "int64",
                "float64",
                "bool",
            }:
                return {
                    "type": "h5",
                    "op": op,
                    "key": widget_session_key,
                    "value": {"error": f"Invalid dtype: {obs_dtype}"},
                }

            adata.obs = adata.obs.reindex(
                columns=[*adata.obs.columns.tolist(), obs_key]
            )
            if obs_dtype == "int64":
                adata.obs[obs_key] = adata.obs[obs_key].fillna(0).astype("int64")
            else:
                adata.obs[obs_key] = adata.obs[obs_key].astype(str(obs_dtype) if obs_dtype is not None else "category")  # type: ignore

            created_for_key = obs_key

        mutated_for_key = None
        if "obs_value" in msg and "lasso_points" in msg and "obsm_key" in msg:
            mutate_obs_by_lasso(
                adata,
                msg["obsm_key"],
                obs_key,
                msg["obs_value"],
                msg["lasso_points"],
                msg.get("filters"),
            )
            mutated_for_key = obs_key

        if "old_obs_value" in msg and "new_obs_value" in msg:
            mutate_obs_by_value(
                adata, obs_key, msg["old_obs_value"], msg["new_obs_value"]
            )
            mutated_for_key = obs_key

        obs, (unique_obs, counts), nrof_obs = get_obs(
            obj_id, adata, obs_key, max_visualization_cells
        )

        return {
            "type": "h5",
            "op": "get_obs",
            "data_type": "h5ad",
            "key": widget_session_key,
            "value": {
                "data": {
                    "fetched_for_key": msg["obs_key"],
                    "mutated_for_key": mutated_for_key,
                    "created_for_key": created_for_key,
                    "values": obs.tolist(),
                    "unique_values": unique_obs.tolist(),
                    "counts": counts.tolist(),
                    "nrof_values": nrof_obs,
                    "dtype": (
                        str(adata.obs[str(obs_key)].dtype)
                        if obs_key in adata.obs
                        else None
                    ),
                },
            },
        }

    if op == "drop_obs":
        if "obs_key" not in msg:
            return {
                "type": "h5",
                "op": op,
                "data_type": "h5ad",
                "key": widget_session_key,
                "value": {"error": "`obs_key` key missing from message"},
            }

        obs_key = msg["obs_key"]
        if obs_key not in adata.obs:
            return {
                "type": "h5",
                "op": op,
                "data_type": "h5ad",
                "key": widget_session_key,
                "value": {"error": "Observation key not found"},
            }

        adata.obs = adata.obs.drop(columns=[obs_key])

        return {
            "type": "h5",
            "op": op,
            "data_type": "h5ad",
            "key": widget_session_key,
            "value": {
                "data": {
                    "dropped_key": obs_key,
                }
            },
        }

    if op == "rename_obs":
        if "old_obs_key" not in msg or "new_obs_key" not in msg:
            return {
                "type": "h5",
                "op": op,
                "data_type": "h5ad",
                "key": widget_session_key,
                "value": {
                    "error": "`old_obs_key` and `new_obs_key` keys missing from message"
                },
            }

        old_obs_key = msg["old_obs_key"]
        new_obs_key = msg["new_obs_key"]

        if old_obs_key not in adata.obs:
            return {
                "type": "h5",
                "op": op,
                "data_type": "h5ad",
                "key": widget_session_key,
                "value": {"error": "Observation key not found"},
            }

        adata.obs = adata.obs.rename(columns={old_obs_key: new_obs_key})

        return {
            "type": "h5",
            "op": op,
            "data_type": "h5ad",
            "key": widget_session_key,
            "value": {
                "data": {
                    "old_key": old_obs_key,
                    "new_key": new_obs_key,
                }
            },
        }

    if op == "fetch_and_process_image":
        if "s3_presigned_url" not in msg or "node_id" not in msg:
            return {
                "type": "h5",
                "op": op,
                "data_type": "h5ad",
                "key": widget_session_key,
                "value": {
                    "error": "`s3_presigned_url` or `node_id` key missing from message"
                },
            }

        image_uri = await fetch_and_process_image(
            msg["node_id"], msg["s3_presigned_url"]
        )
        return {
            "type": "h5",
            "op": op,
            "data_type": "h5ad",
            "key": widget_session_key,
            "value": {
                "data": {"image": image_uri, "fetched_for_node_id": msg.get("node_id")}
            },
        }

    if op == "align_image":
        alignment_is_running = True
        try:
            try:
                image_bytes = pil_image_cache[msg["node_id"]]
            except KeyError:
                return {
                    "type": "h5",
                    "op": op,
                    "data_type": "h5ad",
                    "key": widget_session_key,
                    "value": {
                        "data": {
                            "stage": "fetch_image",
                            "error": (
                                f"attempting to align image from an unprocessed node (nid: {msg['node_id']})"
                            ),
                        }
                    },
                }

            await align_image(
                msg["scatter_data_key"],
                msg["new_scatter_data_key"],
                msg["points_I"],
                msg["points_J"],
                msg["alignment_method"],
                image_bytes,
                adata,
                widget_session_key,
                send,
            )
            alignment_is_running = False
            return {
                "type": "h5",
                "op": op,
                "data_type": "h5ad",
                "key": widget_session_key,
                "value": {"data": {"aligned_obsm_key": msg["new_scatter_data_key"]}},
            }
        finally:
            alignment_is_running = False

    if op == "store_views":
        if "views" not in msg:
            return {
                "type": "h5",
                "op": op,
                "data_type": "h5ad",
                "key": widget_session_key,
                "value": {"error": "`views` key missing from message"},
            }

        adata.uns["latch_views"] = msg["views"]

        return {
            "type": "h5",
            "op": op,
            "data_type": "h5ad",
            "key": widget_session_key,
            "value": {
                "data": {
                    "stored_views": adata.uns["latch_views"],
                },
            },
        }

    raise ValueError(f"Invalid operation: {op}")
