import base64
from io import BytesIO
from typing import Any

import aiohttp
import numpy as np
import pandas as pd
from matplotlib.path import Path
from numpy.typing import NDArray
from PIL import Image

from lplots.h5.utils import auto_install

Image.MAX_IMAGE_PIXELS = None  # TIFFs can be huge, and PIL detects a decompression bomb on many ~400mb examples otherwise

ad = auto_install.ad

ann_data_index_cache: dict[str, tuple[int, NDArray[np.bool_], NDArray[np.int64]]] = {}
ann_data_var_index_cache: dict[str, tuple[NDArray[np.str_], NDArray[np.str_] | None]] = {}


rng = np.random.default_rng()

pil_image_cache: dict[str, bytes] = {}


async def fetch_and_process_image(
    node_id: str,
    s3_presigned_url: str,
    max_width: int = 1536,
    max_height: int = 1536,
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
        img.save(output_buffer, format="PNG", optimize=True)
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
    max_visualization_cells: int = 100000,
) -> tuple[NDArray[np.float32], NDArray[np.str_], bool] | tuple[None, None, bool]:
    if obsm_key not in adata.obsm:
        return None, None, False

    mask = generate_filter_mask(adata, filters if filters is not None else [])
    filtered_n_cells = np.sum(mask)
    recomputed_index = False

    # todo(aidan): intelligent downsampling to preserve outliers / information in general
    if filtered_n_cells > max_visualization_cells:
        filtered_indices = np.where(mask)[0]
        stored_for_max_cells, stored_mask, stored_idxs = ann_data_index_cache.get(obj_id, (None, None, None))
        if stored_mask is None or stored_idxs is None or not np.array_equal(stored_mask, mask) or stored_for_max_cells != max_visualization_cells:
            idxs = rng.choice(filtered_indices, size=max_visualization_cells, replace=False)
            ann_data_index_cache[obj_id] = (max_visualization_cells, mask, idxs)
            recomputed_index = True
        else:
            idxs = stored_idxs
    else:
        stored_for_max_cells, stored_mask, stored_idxs = ann_data_index_cache.get(obj_id, (None, None, None))
        if stored_mask is None or stored_idxs is None or not np.array_equal(stored_mask, mask) or stored_for_max_cells != max_visualization_cells:
            idxs = np.where(mask)[0]
            ann_data_index_cache[obj_id] = (max_visualization_cells, mask, idxs)
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
    max_visualization_cells: int = 100000,
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

    _, _, idxs = ann_data_index_cache[obj_id]

    return obs[idxs], (truncated_unique_obs, counts), len(unique_obs)


def get_obs_vector(
    obj_id: str,
    adata: ad.AnnData,
    var_index: str,
) -> NDArray[np.int64]:
    _, _, idxs = ann_data_index_cache[obj_id]

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
