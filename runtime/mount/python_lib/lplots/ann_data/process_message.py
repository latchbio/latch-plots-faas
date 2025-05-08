import base64
from io import BytesIO
from typing import Any

import aiohttp
import numpy as np
import pandas as pd
from matplotlib.path import Path
from numpy.typing import NDArray
from PIL import Image

from lplots.ann_data import align_image, auto_install

from .. import _inject

Image.MAX_IMAGE_PIXELS = None  # TIFFs can be huge, and PIL detects a decompression bomb on many ~400mb examples otherwise

ad = auto_install.ad

ann_data_index_cache: dict[str, tuple[NDArray[np.bool_], NDArray[np.int64]]] = {}
ann_data_var_index_cache: dict[str, tuple[NDArray[np.str_], NDArray[np.str_] | None]] = {}


max_visualization_cells = 100000
rng = np.random.default_rng()

pil_image_cache: dict[str, bytes] = {}


async def fetch_and_process_image(
    node_id: str,
    s3_presigned_url: str,
    max_width: int = 512,
    max_height: int = 512,
) -> str:
    if node_id in pil_image_cache:
        data = pil_image_cache[node_id]
    else:
        async with aiohttp.ClientSession() as session, session.get(s3_presigned_url) as response:
            if response.status != 200:
                raise ValueError(f"Failed to fetch image. Status: {response.status}")
            data = await response.read()
            pil_image_cache[node_id] = data

    image_data = BytesIO(data)
    with Image.open(image_data) as img:
        img.thumbnail((max_width, max_height))  # `thumbnail` maintains aspect ratio, `resize` does not

        output_buffer = BytesIO()
        img.save(output_buffer, format="PNG")
        output_buffer.seek(0)

        base64_str = base64.b64encode(output_buffer.read()).decode("utf-8")

    return f"data:image/png;base64,{base64_str}"


def generate_filter_mask(
    adata: ad.AnnData,
    filters: list[dict[str, Any]],
) -> NDArray[np.bool_]:
    mask = np.ones(adata.n_obs, dtype=bool)

    for f in filters:
        if f["type"] == "obs":
            key = f["key"]
            if key not in adata.obs:
                continue

            op = f["operation"]
            value = op["value"]
            obs_values = adata.obs[key]

            if op["type"] == "neq":
                if value is None:
                    mask &= ~pd.isna(obs_values)
                else:
                    mask &= obs_values != value
            elif op["type"] == "geq":
                mask &= obs_values >= value
            elif op["type"] == "leq":
                mask &= obs_values <= value
            elif op["type"] == "g":
                mask &= obs_values > value
            elif op["type"] == "l":
                mask &= obs_values < value

        elif f["type"] == "var":
            key = f["key"]
            if key not in adata.var_names:
                continue

            op = f["operation"]
            value = op["value"]

            # todo(aidan): likely too slow?
            var_values = np.asarray(adata[:, key].to_df().iloc[:, 0:].values.ravel())  # noqa: PD011

            if op["type"] == "geq":
                mask &= var_values >= value
            elif op["type"] == "leq":
                mask &= var_values <= value
            elif op["type"] == "g":
                mask &= var_values > value
            elif op["type"] == "l":
                mask &= var_values < value

    return mask


def get_obsm(
    obj_id: str,
    adata: ad.AnnData,
    obsm_key: str,
    filters: list[dict[str, Any]] | None = None,
) -> tuple[NDArray[np.float32], NDArray[np.str_], bool] | tuple[None, None, bool]:
    if obsm_key not in adata.obsm:
        return None, None, False

    mask = generate_filter_mask(adata, filters if filters is not None else [])
    filtered_n_cells = np.sum(mask)
    recomputed_index = False

    # todo(aidan): intelligent downsampling to preserve outliers / information in general
    if filtered_n_cells > max_visualization_cells:
        filtered_indices = np.where(mask)[0]
        stored_mask, stored_idxs = ann_data_index_cache.get(obj_id, (None, None))
        if stored_mask is None or stored_idxs is None or not np.array_equal(stored_mask, mask):
            idxs = rng.choice(filtered_indices, size=max_visualization_cells, replace=False)
            ann_data_index_cache[obj_id] = (mask, idxs)
            recomputed_index = True
        else:
            idxs = stored_idxs
    else:
        stored_mask, stored_idxs = ann_data_index_cache.get(obj_id, (None, None))
        if stored_mask is None or stored_idxs is None or not np.array_equal(stored_mask, mask):
            idxs = np.where(mask)[0]
            ann_data_index_cache[obj_id] = (mask, idxs)
            recomputed_index = True
        else:
            idxs = stored_idxs

    # todo(aidan): allow specifying dimensions to fetch (PCA components is an example where this is useful)
    # todo(aidan): support sparse data?
    obsm = np.asarray(adata.obsm[obsm_key][idxs, :2], dtype=np.float32)  # type: ignore
    index = np.asarray(adata.obs_names[idxs])

    return obsm, index, recomputed_index


def get_obs(
    obj_id: str,
    adata: ad.AnnData,
    obs_key: str,
) -> tuple[NDArray[np.str_], tuple[NDArray[np.str_], NDArray[np.int64]], int]:
    obs = np.asarray(adata.obs[obs_key])

    value_counts = pd.Series(obs).value_counts(dropna=False)
    unique_obs = value_counts.index.to_numpy()
    counts = value_counts.values.astype(np.int64)  # noqa: PD011

    if len(unique_obs) > max_visualization_cells:
        sorted_indices = np.argsort(counts)[::-1][:max_visualization_cells]
        truncated_unique_obs = unique_obs[sorted_indices]
        counts = counts[sorted_indices]
    else:
        truncated_unique_obs = unique_obs

    _, idxs = ann_data_index_cache[obj_id]

    return obs[idxs], (truncated_unique_obs, counts), len(unique_obs)


def get_obs_vector(
    obj_id: str,
    adata: ad.AnnData,
    var_index: str,
) -> NDArray[np.int64]:
    _, idxs = ann_data_index_cache[obj_id]

    # note(aidan): this is ~3+ times faster than using `.to_numpy()` in place of `.values`
    return np.asarray(adata[idxs, var_index].to_df().iloc[:, 0:].values.ravel())  # noqa: PD011


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


def adapt_value_for_dtype(
    value: str | float | int | bool | None,  # noqa: FBT001
    dtype: pd.api.extensions.ExtensionDtype | np.dtype,
) -> str | float | bool | None:
    if value is None:
        return value

    # note(aidan): bools are also numeric dtypes so this needs to happen before the numeric check
    if pd.api.types.is_bool_dtype(dtype):
        if isinstance(value, str):
            return value.lower().strip() in {"true", "1", "yes"}

        return bool(value)

    if pd.api.types.is_numeric_dtype(dtype):
        if "int" in str(dtype).lower():
            return int(value)

        return float(value)

    if pd.api.types.is_datetime64_any_dtype(dtype):
        return str(pd.to_datetime(value, errors="coerce"))

    return str(value)


def mutate_obs_by_lasso(
    adata: ad.AnnData,
    obsm_key: str,
    obs_key: str,
    obs_value: str | float | int | bool | None,  # noqa: FBT001
    lasso_points: list[tuple[int, int]],
    filters: list[dict[str, Any]] | None = None,
) -> None:
    embedding = adata.obsm[obsm_key][:, :2]  # type: ignore

    if len(lasso_points) < 3:
        return

    polygon = Path(lasso_points)
    mask = polygon.contains_points(embedding) & generate_filter_mask(adata, filters if filters is not None else [])

    coerced_obs_value = adapt_value_for_dtype(obs_value, adata.obs[obs_key].dtype)

    if isinstance(adata.obs[obs_key].dtype, pd.CategoricalDtype) and coerced_obs_value not in adata.obs[obs_key].cat.categories and coerced_obs_value is not None:
        adata.obs[obs_key] = adata.obs[obs_key].cat.add_categories(coerced_obs_value)
    adata.obs.loc[mask, obs_key] = coerced_obs_value


def mutate_obs_by_value(
    adata: ad.AnnData,
    obs_key: str,
    old_obs_value: str | float | int | bool | None,  # noqa: FBT001
    new_obs_value: str | float | int | bool | None,  # noqa: FBT001
) -> None:
    coerced_new_obs_value = adapt_value_for_dtype(new_obs_value, adata.obs[obs_key].dtype)

    if isinstance(adata.obs[obs_key].dtype, pd.CategoricalDtype) and coerced_new_obs_value not in adata.obs[obs_key].cat.categories and coerced_new_obs_value is not None:
        adata.obs[obs_key] = adata.obs[obs_key].cat.add_categories(coerced_new_obs_value)

    coerced_old_obs_value = adapt_value_for_dtype(old_obs_value, adata.obs[obs_key].dtype)

    if coerced_old_obs_value is None:
        mask = adata.obs[obs_key].isna()
    else:
        mask = adata.obs[obs_key] == coerced_old_obs_value

    if coerced_new_obs_value is None:
        adata.obs.loc[mask, obs_key] = None
    else:
        adata.obs.loc[mask, obs_key] = coerced_new_obs_value

    if isinstance(adata.obs[obs_key].dtype, pd.CategoricalDtype) and coerced_old_obs_value is not None and coerced_old_obs_value in adata.obs[obs_key].cat.categories:
        adata.obs[obs_key] = adata.obs[obs_key].cat.remove_categories(coerced_old_obs_value)


async def handle_ann_data_widget_message(
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

    widget_session_key = msg["key"]
    widget_state: dict[str, Any] = msg["state"]

    if "obj_id" not in widget_state:
        return {
            "type": "ann_data",
            "key": widget_session_key,
            "value": {
                "error": "Widget state does not contain obj_id",
            },
        }

    obj_id = widget_state["obj_id"]
    if obj_id not in _inject.kernel.ann_data_objects:
        return {
            "type": "ann_data",
            "key": widget_session_key,
            "value": {
                "error": f"AnnData object with ID {obj_id} not found in cache",
            },
        }

    adata: ad.AnnData = _inject.kernel.ann_data_objects[obj_id]

    if "op" not in msg or msg["op"] not in {"init_data", "get_obsm_options", "get_obsm", "get_obs_options", "get_obs", "get_counts_column", "mutate_obs", "drop_obs", "rename_obs", "fetch_and_process_image"}:
        return {
            "type": "ann_data",
            "key": widget_session_key,
            "value": {
                "error": f"Invalid operation: {msg.get('op', '`op` key missing from message')}",
            },
        }

    op = msg["op"]

    if op == "init_data":
        init_obsm_key = msg.get("obsm_key")
        possible_obsm_keys = adata.obsm_keys()
        if init_obsm_key is None:
            for key in possible_obsm_keys:
                if "umap" in key.lower():
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
            obsm, index, recomputed_index = get_obsm(obj_id, adata, init_obsm_key, filters)

        obs = None
        unique_obs = None
        nrof_obs = None
        counts = None
        if init_obs_key is not None and init_obs_key in adata.obs:
            obs, (unique_obs, counts), nrof_obs = get_obs(obj_id, adata, init_obs_key)

        gene_column = None
        if init_var_key is not None and init_obs_key is None and init_var_key in adata.var_names:
            gene_column = get_obs_vector(obj_id, adata, init_var_key)

        var_index, var_names = get_var_index(obj_id, adata)

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_session_key,
            "value": {
                "data": {
                    # display info
                    "num_obs": adata.n_obs,
                    "num_vars": adata.n_vars,

                    # options
                    "possible_obs_keys": possible_obs_keys,
                    "possible_obs_keys_types": [str(adata.obs[key].dtype) for key in possible_obs_keys],
                    "possible_obsm_keys": possible_obsm_keys,

                    # init state with these
                    "init_obs_key": init_obs_key,
                    "init_obsm_key": init_obsm_key,
                    "init_recomputed_index": recomputed_index,
                    "init_obsm_values": obsm.tolist() if obsm is not None else None,
                    "init_obsm_index": index.tolist() if index is not None else None,
                    "init_obsm_filters": filters,
                    "init_obs_values": obs.tolist() if obs is not None else None,
                    "init_obs_unique_values": unique_obs.tolist() if unique_obs is not None else None,
                    "init_obs_counts": counts.tolist() if counts is not None else None,
                    "init_obs_nrof_values": nrof_obs,

                    # var info
                    "init_var_index": var_index.tolist(),
                    "init_var_names": var_names.tolist() if var_names is not None else None,

                    # var color by info
                    "init_var_values": gene_column.tolist() if gene_column is not None else None,
                    "init_var_key": init_var_key if init_var_key is not None else None,
                }
            },
        }

    if op == "get_obsm_options":
        return {
            "type": "ann_data",
            "op": op,
            "key": widget_session_key,
            "value": {"data": list(adata.obsm.keys())},
        }

    if op == "get_obsm":
        if "obsm_key" not in msg or msg["obsm_key"] not in adata.obsm:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_session_key,
                "value": {
                    "error": "Obsm not found" if "obsm_key" in msg else "`obsm_key` key missing from message"
                },
            }

        filters = msg.get("filters")
        obsm, index, recomputed_index = get_obsm(obj_id, adata, msg["obsm_key"], filters)
        if obsm is None or index is None:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_session_key,
                "value": {"error": "Failed to get obsm"},
            }

        obs = None
        fetched_for_obs_key = None
        unique_obs = None
        counts = None
        nrof_obs = None
        gene_column = None
        fetched_for_var_key = None
        if "colored_by_type" in msg and "colored_by_key" in msg:
            if msg["colored_by_type"] == "obs" and msg["colored_by_key"] in adata.obs:
                obs, (unique_obs, counts), nrof_obs = get_obs(obj_id, adata, msg["colored_by_key"])
                fetched_for_obs_key = msg["colored_by_key"]
            elif msg["colored_by_type"] == "var" and msg["colored_by_key"] in adata.var_names:
                gene_column = get_obs_vector(obj_id, adata, msg["colored_by_key"])
                fetched_for_var_key = msg["colored_by_key"]

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_session_key,
            "value": {
                "data": {
                    "fetched_for_key": msg["obsm_key"],
                    "obsm": obsm.tolist(),
                    "index": index.tolist(),
                    "recomputed_index": recomputed_index,
                    "fetched_for_obs_key": fetched_for_obs_key,
                    "fetched_for_var_key": fetched_for_var_key,
                    "fetched_for_filters": filters,
                    "values": obs.tolist() if obs is not None else None,
                    "unique_values": unique_obs.tolist() if unique_obs is not None else None,
                    "counts": counts.tolist() if counts is not None else None,
                    "nrof_values": nrof_obs,
                    "var_values": gene_column.tolist() if gene_column is not None else None,
                },
            },
        }

    if op == "get_obs_options":
        return {
            "type": "ann_data",
            "op": op,
            "key": widget_session_key,
            "value": {
                "data": list(adata.obs.keys()),
            },
        }

    if op == "get_obs":
        if "obs_key" not in msg or msg["obs_key"] not in adata.obs:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_session_key,
                "value": {
                    "error": "Observation not found" if "obs_key" in msg else "`obs_key` key missing from message"
                },
            }

        obs, (unique_obs, counts), nrof_obs = get_obs(obj_id, adata, msg["obs_key"])

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_session_key,
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
                "key": widget_session_key,
                "value": {
                    "error": "Variable not found" if "var_index" in msg else "`var_index` key missing from message"
                },
            }

        gene_column = get_obs_vector(obj_id, adata, msg["var_index"])

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_session_key,
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
                "key": widget_session_key,
                "value": {"error": "`obs_key` key missing from message"},
            }

        obs_key = str(msg["obs_key"])
        created_for_key = None
        if obs_key not in adata.obs:
            obs_dtype = msg.get("obs_dtype")
            if obs_dtype is not None and obs_dtype not in {"category", "int64", "float64", "bool"}:
                return {
                    "type": "ann_data",
                    "op": op,
                    "key": widget_session_key,
                    "value": {"error": f"Invalid dtype: {obs_dtype}"},
                }

            adata.obs = adata.obs.reindex(columns=[*adata.obs.columns.tolist(), obs_key])
            if obs_dtype == "int64":
                adata.obs[obs_key] = adata.obs[obs_key].fillna(0).astype("int64")
            else:
                adata.obs[obs_key] = adata.obs[obs_key].astype(str(obs_dtype) if obs_dtype is not None else "category")  # type: ignore

            created_for_key = obs_key

        mutated_for_key = None
        if "obs_value" in msg and "lasso_points" in msg and "obsm_key" in msg:
            mutate_obs_by_lasso(adata, msg["obsm_key"], obs_key, msg["obs_value"], msg["lasso_points"], msg.get("filters"))
            mutated_for_key = obs_key

        if "old_obs_value" in msg and "new_obs_value" in msg:
            mutate_obs_by_value(adata, obs_key, msg["old_obs_value"], msg["new_obs_value"])
            mutated_for_key = obs_key

        obs, (unique_obs, counts), nrof_obs = get_obs(obj_id, adata, obs_key)

        return {
            "type": "ann_data",
            "op": "get_obs",
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
                    "dtype": str(adata.obs[str(obs_key)].dtype) if obs_key in adata.obs else None,
                },
            },
        }

    if op == "drop_obs":
        if "obs_key" not in msg:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_session_key,
                "value": {"error": "`obs_key` key missing from message"},
            }

        obs_key = msg["obs_key"]
        if obs_key not in adata.obs:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_session_key,
                "value": {"error": "Observation key not found"},
            }

        adata.obs = adata.obs.drop(columns=[obs_key])

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_session_key,
            "value": {"data": {
                "dropped_key": obs_key,
            }},
        }

    if op == "rename_obs":
        if "old_obs_key" not in msg or "new_obs_key" not in msg:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_session_key,
                "value": {"error": "`old_obs_key` and `new_obs_key` keys missing from message"},
            }

        old_obs_key = msg["old_obs_key"]
        new_obs_key = msg["new_obs_key"]

        if old_obs_key not in adata.obs:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_session_key,
                "value": {"error": "Observation key not found"},
            }

        adata.obs = adata.obs.rename(columns={old_obs_key: new_obs_key})

        return {
            "type": "ann_data",
            "op": op,
            "key": widget_session_key,
            "value": {"data": {
                "old_key": old_obs_key,
                "new_key": new_obs_key,
            }},
        }

    if op == "fetch_and_process_image":
        if "s3_presigned_url" not in msg or "node_id" not in msg:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_session_key,
                "value": {"error": "`s3_presigned_url` or `node_id` key missing from message"},
            }

        image_uri = await fetch_and_process_image(msg["node_id"], msg["s3_presigned_url"])
        return {
            "type": "ann_data",
            "op": op,
            "key": widget_session_key,
            "value": {"data": {"image": image_uri, "fetched_for_node_id": msg.get("node_id")}},
        }

    if op == "align_image":
        try:
            image_bytes = pil_image_cache[msg["node_id"]]
        except KeyError:
            return {
                "type": "ann_data",
                "op": op,
                "key": widget_session_key,
                "value": {"error": {f"attempting to align image from an unprocessed node (nid: {msg['node_id']})"}},
            }
        aligned_obs_key = align_image(
                msg["points_I"],
                msg["points_J"],
                msg["alignment_method"],
                image_bytes
        )
        return {
            "type": "ann_data",
            "op": op,
            "key": aligned_obs_key,
            "value": {"data": {"aligned_obs_key": aligned_obs_key}},
        }

    raise ValueError(f"Invalid operation: {op}")
