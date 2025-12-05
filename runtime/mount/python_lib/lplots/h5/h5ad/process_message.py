import asyncio
from collections.abc import Awaitable, Callable
from typing import Any, overload

import numpy as np
from numpy.typing import NDArray

from ..utils import auto_install
from ..utils.align import align_image
from .ops import (
    Context,
    fetch_and_process_image,
    get_var_index,
    mutate_obs_by_lasso,
    mutate_obs_by_value,
    pil_image_cache,
)

ad = auto_install.ad

alignment_is_running = False

contexts: dict[str, Context] = {}


async def process_h5ad_request(
    msg: dict[str, Any],
    *,
    widget_session_key: str,
    obj_id: str,
    send: Callable[[object], Awaitable[None]],
) -> dict[str, Any] | None:
    global alignment_is_running

    op = msg.get("op")
    if op is None:
        return {
            "type": "h5",
            "key": widget_session_key,
            "value": {"error": ("Invalid operation: `op` key missing from message")},
        }

    ctx = contexts.get(obj_id)
    if ctx is None:
        ctx = Context(id=obj_id)
        contexts[obj_id] = ctx
    adata = ctx.adata

    @overload
    def make_response(*, data: Any) -> dict[str, Any]: ...
    @overload
    def make_response(*, error: str) -> dict[str, Any]: ...
    def make_response(
        *, data: Any | None = None, error: str | None = None
    ) -> dict[str, Any]:
        value = {"data": data} if data is not None else {"error": error}
        return {
            "type": "h5",
            "op": op,
            "data_type": "h5ad",
            "key": widget_session_key,
            "value": value,
        }

    max_cells = msg.get("max_visualization_cells", 100000)
    assert isinstance(max_cells, int)

    match op:
        case "init_data":
            init_obsm_key = msg.get("obsm_key")
            possible_obsm_keys = adata.obsm_keys()
            if init_obsm_key is None:
                for key in possible_obsm_keys:
                    if "umap" in key.lower() or "spatial" in key.lower():
                        init_obsm_key = key
                        break

            if (
                init_obsm_key is None or init_obsm_key not in possible_obsm_keys
            ) and len(possible_obsm_keys) > 0:
                init_obsm_key = possible_obsm_keys[0]

            possible_obs_keys = adata.obs_keys()

            init_obs_key = msg.get("obs_key")
            assert isinstance(init_obs_key, str | None)

            init_var_key = msg.get("var_key")
            assert isinstance(init_var_key, str | None)

            if init_obs_key is None and init_var_key is None:
                for key in possible_obs_keys:
                    if "cell" in key.lower() and "type" in key.lower():
                        init_obs_key = key
                        break

            if (
                init_obs_key is None
                and init_var_key is None
                and len(possible_obs_keys) > 0
            ):
                init_obs_key = possible_obs_keys[0]

            filters = msg.get("filters")
            recomputed_index = ctx.compute_index(filters=filters, max_cells=max_cells)

            obsm = None
            if init_obsm_key is not None:
                obsm = ctx.get_obsm(init_obsm_key)

            obs = None
            if init_obs_key is not None:
                obs = ctx.get_obs(init_obs_key, max_cells=max_cells)

            gene_column = None
            if init_var_key is not None and init_obs_key is None:
                gene_column = ctx.get_obs_vector(init_var_key)

            var_index, var_names = get_var_index(obj_id, adata)

            global alignment_is_running

            return make_response(
                data={
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
                    "init_obsm_values": obsm.data.tolist()
                    if obsm is not None
                    else None,
                    "init_obsm_index": obsm.index.tolist()
                    if obsm is not None
                    else None,
                    "init_obsm_filters": filters,
                    "init_obs_values": obs.data.tolist() if obs is not None else None,
                    "init_obs_unique_values": (
                        obs.top_values.tolist() if obs is not None else None
                    ),
                    "init_obs_counts": obs.top_value_counts.tolist()
                    if obs is not None
                    else None,
                    "init_obs_nrof_values": obs.total_unique
                    if obs is not None
                    else None,
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
                    # images
                    "init_images": adata.uns.get("latch_images", {}),
                }
            )

        case "get_obsm_options":
            return make_response(data=list(adata.obsm.keys()))

        case "get_obsm":
            if "obsm_key" not in msg:
                return make_response(error="`obsm_key` key missing from message")

            filters = msg.get("filters")
            recomputed_index = ctx.compute_index(filters=filters, max_cells=max_cells)

            obsm = ctx.get_obsm(msg["obsm_key"])
            if obsm is None:
                return make_response(error="Obsm not found")

            obs = None
            fetched_for_obs_key = None
            gene_columns: list[NDArray[np.int64]] | None = None
            fetched_for_var_keys: list[str] | None = None
            if "colored_by_type" in msg and "colored_by_key" in msg:
                if msg["colored_by_type"] == "obs":
                    obs = ctx.get_obs(msg["colored_by_key"], max_cells=max_cells)

                    if obs is not None:
                        fetched_for_obs_key = msg["colored_by_key"]
                elif msg["colored_by_type"] == "var":
                    fetched_for_var_keys = [
                        x for x in msg["colored_by_key"] if x in adata.var_names
                    ]

                    gene_columns = []
                    for x in fetched_for_var_keys:
                        cur = ctx.get_obs_vector(x)
                        assert cur is not None

                        gene_columns.append(cur)

            return make_response(
                data={
                    "fetched_for_key": msg["obsm_key"],
                    "obsm": obsm.data.tolist(),
                    "index": obsm.index.tolist(),
                    "recomputed_index": recomputed_index,
                    "fetched_for_obs_key": fetched_for_obs_key,
                    "fetched_for_var_keys": fetched_for_var_keys,
                    "fetched_for_filters": filters,
                    "values": obs.data.tolist() if obs is not None else None,
                    "unique_values": (
                        obs.top_values.tolist() if obs is not None else None
                    ),
                    "counts": obs.top_value_counts.tolist()
                    if obs is not None
                    else None,
                    "nrof_values": obs.total_unique if obs is not None else None,
                    "var_values": [col.tolist() for col in gene_columns]
                    if gene_columns is not None
                    else None,
                }
            )

        case "get_obs_options":
            return make_response(data=list(adata.obs.keys()))

        case "get_obs":
            if "obs_key" not in msg:
                return make_response(error="`obs_key` key missing from message")

            obs = ctx.get_obs(msg["obs_key"], max_cells=max_cells)
            if obs is None:
                return make_response(error="Observation not found")

            return make_response(
                data={
                    "fetched_for_key": msg["obs_key"],
                    "values": obs.data.tolist(),
                    "unique_values": obs.top_values.tolist(),
                    "counts": obs.top_value_counts.tolist(),
                    "nrof_values": obs.total_unique,
                }
            )

        case "get_counts_column":
            if "var_index" not in msg:
                return make_response(error="`var_index` key missing from message")

            gene_column = ctx.get_obs_vector(msg["var_index"])
            if gene_column is None:
                return make_response(error="Variable not found")

            return make_response(
                data={
                    "fetched_for_key": msg["var_index"],
                    "values": gene_column.tolist(),
                }
            )

        case "mutate_obs":
            if "obs_key" not in msg:
                return make_response(error="`obs_key` key missing from message")

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
                    return make_response(error=f"Invalid dtype: {obs_dtype}")

                adata.obs = adata.obs.reindex(
                    columns=[*adata.obs.columns.tolist(), obs_key]
                )
                if obs_dtype == "int64":
                    adata.obs[obs_key] = adata.obs[obs_key].fillna(0).astype("int64")
                else:
                    adata.obs[obs_key] = adata.obs[obs_key].astype(
                        str(obs_dtype) if obs_dtype is not None else "category"
                    )  # type: ignore

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

            obs = ctx.get_obs(obs_key, max_cells=max_cells)
            assert obs is not None

            return make_response(
                data={
                    "fetched_for_key": msg["obs_key"],
                    "mutated_for_key": mutated_for_key,
                    "created_for_key": created_for_key,
                    "values": obs.data.tolist(),
                    "unique_values": obs.top_values.tolist(),
                    "counts": obs.top_value_counts.tolist(),
                    "nrof_values": obs.total_unique,
                    "dtype": (
                        str(adata.obs[str(obs_key)].dtype)
                        if obs_key in adata.obs
                        else None
                    ),
                }
            )

        case "drop_obs":
            if "obs_key" not in msg:
                return make_response(error="`obs_key` key missing from message")

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

            return make_response(data={"dropped_key": obs_key})

        case "rename_obs":
            if "old_obs_key" not in msg or "new_obs_key" not in msg:
                return make_response(
                    error="`old_obs_key` and `new_obs_key` keys missing from message"
                )

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

            return make_response(data={"old_key": old_obs_key, "new_key": new_obs_key})

        case "fetch_and_process_image":
            if (
                "s3_presigned_url" not in msg
                or "node_id" not in msg
                or "id" not in msg
                or "name" not in msg
                or "transformations" not in msg
            ):
                return make_response(
                    error="`s3_presigned_url` or `node_id` key missing from message"
                )

            image_uri = await fetch_and_process_image(
                msg["node_id"],
                msg["s3_presigned_url"],
                msg.get("max_width", 1536),
                msg.get("max_height", 1536),
            )

            image_id = msg["id"]

            images = adata.uns.get("latch_images", {})
            image_data = images.get(image_id, {})

            image_data["b64_image"] = image_uri
            image_data["node_id"] = msg["node_id"]
            image_data["id"] = image_id
            image_data["name"] = msg["name"]
            image_data["transformations"] = msg["transformations"]
            images[image_id] = image_data
            adata.uns["latch_images"] = images

            return make_response(
                data={"image_data": image_data, "fetched_for_image_id": image_id}
            )

        case "align_image":
            try:
                image_bytes = pil_image_cache[msg["node_id"]]
            except KeyError:
                return make_response(
                    data={
                        "stage": "fetch_image",
                        "error": (
                            f"attempting to align image from an unprocessed node (nid: {msg['node_id']})"
                        ),
                    }
                )

            if alignment_is_running:
                return make_response(
                    data={"stage": "validation", "error": "alignment already running"}
                )

            alignment_is_running = True

            async def run_alignment() -> None:
                global alignment_is_running

                try:
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
                finally:
                    alignment_is_running = False

            # todo(aidan): our websocket handler processes messages in serial
            asyncio.create_task(run_alignment())

            return None

        case "store_image_transformation":
            if "image_transformation" not in msg or "id" not in msg:
                return make_response(
                    error="`image_transformation` or `id` key missing from message"
                )

            image_data = adata.uns.get("latch_images", {}).get(msg["id"])

            if image_data is None:
                return make_response(error="Image data not found")

            image_data["transformations"] = msg["image_transformation"]
            adata.uns["latch_images"][msg["id"]] = image_data

            return make_response(
                data={"updated_image_transformations_for_id": msg["id"]}
            )

        case "remove_image":
            if "id" not in msg:
                return make_response(error="`id` key missing from message")

            adata.uns["latch_images"].pop(msg["id"], None)

            return make_response(data={"removed_image_id": msg["id"]})

        case "store_views":
            if "views" not in msg:
                return make_response(error="`views` key missing from message")

            adata.uns["latch_views"] = msg["views"]

            return make_response(data={"stored_views": adata.uns["latch_views"]})

        case _:
            return make_response(error=f"Invalid operation: {op}")
