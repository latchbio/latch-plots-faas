from dataclasses import dataclass
from typing import Any, Literal, NotRequired, TypedDict

import pandas as pd
from latch.ldata.path import LPath
from latch.registry.table import Table

from .. import _inject
from ..reactive import Signal
from . import _emit, _state
from .shared import FormInputAppearance

DataSourceType = Literal["ldata", "dataframe", "viewer", "registry"]


class LDataDataSource(TypedDict):
    node_id: str
    type: Literal["ldata"]


class DataFrameDataSource(TypedDict):
    key: str
    type: Literal["dataframe"]


class RegistryDataSource(TypedDict):
    table_id: str
    type: Literal["registry"]


class ValueViewerDataSource(TypedDict):
    viewer_id: str
    type: Literal["viewer"]


DataSourceValue = (
    LDataDataSource | DataFrameDataSource | RegistryDataSource | ValueViewerDataSource
)


class DataSourceSelectState(
    _emit.WidgetState[Literal["tabular_datasource_select"], DataSourceValue]
):
    label: str
    readonly: bool
    default: NotRequired[DataSourceValue | None]
    appearance: NotRequired[FormInputAppearance | None]
    required: bool


def get_ldata_df(node_id: Any) -> pd.DataFrame | None:
    if node_id is None or not isinstance(node_id, str):
        return None

    lpath = LPath(f"latch://{node_id}.node")

    name = lpath.name()
    if name is None:
        return None

    if name.endswith(".csv"):
        return pd.read_csv(lpath.download())
    if name.endswith(".xlsx"):
        return pd.read_excel(lpath.download())
    if name.endswith(".tsv"):
        return pd.read_csv(lpath.download(), sep="\t")
    return None


def get_kernel_df(df_id: Any) -> pd.DataFrame | None:
    if df_id is None:
        return None

    g = _inject.kernel.k_globals.get_signal(df_id)
    if g is None:
        return None

    return g()


def get_registry_df(table_id: Any) -> pd.DataFrame | None:
    if table_id is None:
        return None

    return Table(id=table_id).get_dataframe()


@dataclass(frozen=True, kw_only=True)
class TabularDatasourcePicker:
    _key: str
    _state: DataSourceSelectState
    _signal: Signal[object]

    def _value(self, val: object) -> pd.DataFrame | None:
        if not isinstance(val, dict) or not all(
            key in val for key in DataSourceValue.__annotations__
        ):
            val = self._state.get("default")
            # todo(manske): validate default
            if val is None:
                return None

        res_type = val.get("type")
        if res_type == "ldata":
            node_id = val.get("node_id")
            return get_ldata_df(node_id)

        if res_type == "dataframe":
            df_id = val.get("key")
            return get_kernel_df(df_id)

        if res_type == "viewer":
            viewer_id = val.get("viewer_id")
            if viewer_id is None:
                return None

            key, key_type = _inject.kernel.viewer_cell_selections.get(viewer_id)
            if key_type == "ldata_node_id":
                return get_ldata_df(key)

            if key_type == "key":
                return get_kernel_df(key)

            if key_type == "registry_tabe_id":
                return get_registry_df(key)

            if key_type == "url":
                return _inject.kernel.url_dataframes.get(key)

        if res_type == "registry":
            table_id = val.get("table_id")
            return get_registry_df(table_id)

        return None

    @property
    def value(self) -> pd.DataFrame | None:
        val = self._signal()
        return self._value(val)

    def sample(self) -> pd.DataFrame | None:
        val = self._signal.sample()
        return self._value(val)


def w_datasource_picker(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    default: DataSourceValue | None = None,
    appearance: FormInputAppearance | None = None,
    required: bool = False,
) -> TabularDatasourcePicker:
    key = _state.use_state_key(key=key)

    res = TabularDatasourcePicker(
        _key=key,
        _state={
            "type": "tabular_datasource_select",
            "label": label,
            "readonly": readonly,
            "default": default,
            "appearance": appearance,
            "required": required,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
