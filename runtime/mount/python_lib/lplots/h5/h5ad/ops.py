import base64
from dataclasses import dataclass
from io import BytesIO
from typing import Any, Literal, TypedDict

import aiohttp
import numpy as np
import pandas as pd
from matplotlib.path import Path
from numpy.typing import NDArray
from PIL import Image

from ... import _inject
from ..utils import auto_install

Image.MAX_IMAGE_PIXELS = None  # TIFFs can be huge, and PIL detects a decompression bomb on many ~400mb examples otherwise

ad = auto_install.ad


@dataclass(kw_only=True)
class IndexEntry:
    max_cells: int
    mask: NDArray[np.bool_]
    index: NDArray[np.int64]


@dataclass(kw_only=True)
class VarIndexEntry:
    index: NDArray[np.str_]
    names: NDArray[np.str_]


@dataclass(kw_only=True)
class ObsmData:
    data: NDArray[np.float32]
    index: NDArray[np.str_]


@dataclass(kw_only=True)
class ObsData:
    data: NDArray[np.str_ | np.int64]

    top_values: NDArray[np.str_ | np.int64]
    top_value_counts: NDArray[np.int64]

    total_unique: int

    min: np.int64 | None
    max: np.int64 | None


class ColorPalettes(TypedDict):
    continuous: list[str]
    categorical: list[str]


@dataclass(kw_only=True)
class Context:
    id: str

    _index: IndexEntry | None = None
    _var_index: VarIndexEntry | None = None

    @property
    def adata(self) -> ad.AnnData:
        return _inject.kernel.ann_data_objects[self.id]

    @property
    def index(self) -> NDArray[np.int64]:
        assert self._index is not None, (
            f"h5 not initialized correctly: no index calculated (id={self.id!r})"
        )
        return self._index.index

    def compute_index(
        self, *, filters: list[dict[str, Any]] | None = None, max_cells: int
    ) -> bool:
        if filters is None:
            filters = []

        mask = generate_filter_mask(self.adata, filters)

        cached = self._index
        if (
            cached is not None
            and np.array_equal(cached.mask, mask)
            and cached.max_cells != max_cells
        ):
            return False

        visible_cells = np.sum(mask)

        index = np.where(mask)[0]
        if visible_cells > max_cells:
            # todo(aidan): intelligent downsampling to preserve outliers / information in general
            index = rng.choice(index, size=max_cells, replace=False)

        self._index = IndexEntry(max_cells=max_cells, mask=mask, index=index)

        return True

    def get_obsm(self, key: str) -> ObsmData | None:
        if key not in self.adata.obsm:
            return None

        # todo(aidan): allow specifying dimensions to fetch (PCA components is an example where this is useful)
        # todo(aidan): support sparse data?
        data = self.adata.obsm[key]
        if isinstance(data, pd.DataFrame):
            res = data.iloc[self.index, :2].to_numpy(dtype=np.float32, copy=False)
        else:
            res = np.asarray(data[self.index, :2], dtype=np.float32)

        return ObsmData(data=res, index=np.asarray(self.adata.obs_names[self.index]))

    def get_obs(self, key: str, *, max_cells: int) -> ObsData | None:
        if key not in self.adata.obs:
            return None

        obs = np.asarray(self.adata.obs[key])

        value_counts = pd.Series(obs).value_counts(dropna=False)
        unique_obs = value_counts.index.to_numpy()
        counts = value_counts.values.astype(np.int64)  # noqa: PD011

        top = unique_obs
        if len(top) > max_cells:
            # todo(maximsmol): this argsort seems useless cause value_counts already sorts
            sorted_indices = np.argsort(counts)[::-1][:max_cells]

            # todo(maximsmol): I think we should subset the entire series with `.loc` before we decompose it into columns
            top = top[sorted_indices]
            counts = counts[sorted_indices]

        p0 = None
        p100 = None
        if pd.api.types.is_numeric_dtype(obs.dtype):
            p0, p100 = np.quantile(obs, [0, 1])

        return ObsData(
            data=obs[self.index],
            top_values=top,
            top_value_counts=counts,
            total_unique=len(unique_obs),
            min=p0,
            max=p100,
        )

    def get_obs_vector(self, var: str) -> NDArray[np.int64] | None:
        if var not in self.adata.var_names:
            return None

        # note(aidan): this is ~3+ times faster than using `.to_numpy()` in place of `.values`
        return np.asarray(
            self.adata[self.index, var].to_df().iloc[:, 0:].values.ravel()  # noqa: PD011
        )

    def get_vars_range(self, keys: list[str]) -> tuple[int, int] | None:
        data: list[NDArray[np.int64]] = []
        for k in keys:
            if k not in self.adata.var_names:
                continue

            datum = self.adata[:, k].X
            if datum is None:
                continue

            data.append(datum)

        if len(data) == 0:
            return None

        measures = np.hstack(data)
        means = np.mean(measures, axis=1)
        p0, p100 = np.quantile(means, [0, 1])
        return p0, p100

    def export_png(
        self,
        *,
        obsm_key: str,
        data: list[dict[str, Any]],
        layout: object,
        color_palettes: ColorPalettes,
        color_by: tuple[Literal["obs"], str] | tuple[Literal["var"], list[str]] | None,
    ) -> bytes:
        import plotly.graph_objects as go

        df = self.adata.obsm[obsm_key][self.index]
        data[0]["x"] = df[:, 0]
        data[0]["y"] = df[:, 1]

        fig = go.Figure(data=data, layout=layout)
        return fig.to_image(format="png")


rng = np.random.default_rng()

pil_image_cache: dict[str, bytes] = {}


async def fetch_and_process_image(
    node_id: str, s3_presigned_url: str, max_width: int = 1536, max_height: int = 1536
) -> str:
    if node_id in pil_image_cache:
        data = pil_image_cache[node_id]
    else:
        async with (
            aiohttp.ClientSession() as session,
            session.get(s3_presigned_url) as response,
        ):
            if response.status != 200:
                raise ValueError(f"Failed to fetch image. Status: {response.status}")
            data = await response.read()
            pil_image_cache[node_id] = data

    image_data = BytesIO(data)
    with Image.open(image_data) as img:
        img.thumbnail((
            max_width,
            max_height,
        ))  # `thumbnail` maintains aspect ratio, `resize` does not

        output_buffer = BytesIO()
        img.save(output_buffer, format="PNG", optimize=True)
        output_buffer.seek(0)

        base64_str = base64.b64encode(output_buffer.read()).decode("utf-8")

    return f"data:image/png;base64,{base64_str}"


def generate_filter_mask(
    adata: ad.AnnData, filters: list[dict[str, Any]]
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
            keys = f.get("keys", [])
            if not keys:
                # old key format
                key = f.get("key")
                if key:
                    keys = [key]

            valid_keys = [k for k in keys if k in adata.var_names]
            if not valid_keys:
                continue

            op = f["operation"]
            value = op["value"]

            if len(valid_keys) == 1:
                var_values = np.asarray(
                    adata[:, valid_keys[0]].to_df().iloc[:, 0:].values.ravel()
                )  # noqa: PD011
            else:
                gene_values = [
                    np.asarray(adata[:, key].to_df().iloc[:, 0:].values.ravel())  # noqa: PD011
                    for key in valid_keys
                ]
                var_values = np.mean(gene_values, axis=0)

            if op["type"] == "geq":
                mask &= var_values >= value
            elif op["type"] == "leq":
                mask &= var_values <= value
            elif op["type"] == "g":
                mask &= var_values > value
            elif op["type"] == "l":
                mask &= var_values < value

    return mask


ann_data_var_index_cache: dict[
    str, tuple[NDArray[np.str_], NDArray[np.str_] | None]
] = {}


def get_var_index(
    obj_id: str, adata: ad.AnnData
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
    obsm_raw = adata.obsm[obsm_key]
    if isinstance(obsm_raw, pd.DataFrame):
        embedding = obsm_raw.iloc[:, :2].to_numpy()
    else:
        embedding = np.asarray(obsm_raw[:, :2])

    if len(lasso_points) < 3:
        return

    polygon = Path(lasso_points)
    mask = polygon.contains_points(embedding) & generate_filter_mask(
        adata, filters if filters is not None else []
    )

    coerced_obs_value = adapt_value_for_dtype(obs_value, adata.obs[obs_key].dtype)

    if (
        isinstance(adata.obs[obs_key].dtype, pd.CategoricalDtype)
        and coerced_obs_value not in adata.obs[obs_key].cat.categories
        and coerced_obs_value is not None
    ):
        adata.obs[obs_key] = adata.obs[obs_key].cat.add_categories(coerced_obs_value)
    adata.obs.loc[mask, obs_key] = coerced_obs_value


def mutate_obs_by_value(
    adata: ad.AnnData,
    obs_key: str,
    old_obs_value: str | float | int | bool | None,  # noqa: FBT001
    new_obs_value: str | float | int | bool | None,  # noqa: FBT001
) -> None:
    coerced_new_obs_value = adapt_value_for_dtype(
        new_obs_value, adata.obs[obs_key].dtype
    )

    if (
        isinstance(adata.obs[obs_key].dtype, pd.CategoricalDtype)
        and coerced_new_obs_value not in adata.obs[obs_key].cat.categories
        and coerced_new_obs_value is not None
    ):
        adata.obs[obs_key] = adata.obs[obs_key].cat.add_categories(
            coerced_new_obs_value
        )

    coerced_old_obs_value = adapt_value_for_dtype(
        old_obs_value, adata.obs[obs_key].dtype
    )

    if coerced_old_obs_value is None:
        mask = adata.obs[obs_key].isna()
    else:
        mask = adata.obs[obs_key] == coerced_old_obs_value

    if coerced_new_obs_value is None:
        adata.obs.loc[mask, obs_key] = None
    else:
        adata.obs.loc[mask, obs_key] = coerced_new_obs_value

    if (
        isinstance(adata.obs[obs_key].dtype, pd.CategoricalDtype)
        and coerced_old_obs_value is not None
        and coerced_old_obs_value in adata.obs[obs_key].cat.categories
    ):
        adata.obs[obs_key] = adata.obs[obs_key].cat.remove_categories(
            coerced_old_obs_value
        )
