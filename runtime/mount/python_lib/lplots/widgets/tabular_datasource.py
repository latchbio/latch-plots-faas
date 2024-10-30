from dataclasses import dataclass
from typing import Any, Literal, NotRequired, TypedDict

import latch_data_validation
import latch_data_validation.data_validation
import pandas as pd
from latch.ldata.path import LPath
from latch.registry.table import Table

from .. import _inject
from ..reactive import Signal
from . import _emit, _state
from .shared import FormInputAppearance

DataSourceType = Literal["ldata", "dataframe", "registry"]


class LDataDataSource(TypedDict):
    path: str
    type: Literal["ldata"]


class DataFrameDataSource(TypedDict):
    df_id: str
    type: Literal["dataframe"]


class RegistryDataSource(TypedDict):
    table_id: str
    type: Literal["registry"]


DataSourceValue = LDataDataSource | DataFrameDataSource | RegistryDataSource


class DataSourceSelectState(
    _emit.WidgetState[Literal["tabular_datasource_select"], DataSourceValue]
):
    label: str
    readonly: bool
    default: NotRequired[DataSourceValue | None]
    appearance: NotRequired[FormInputAppearance | None]
    required: bool


@dataclass(frozen=True, kw_only=True)
class TabularDatasourcePicker:
    _key: str
    _state: DataSourceSelectState
    _signal: Signal[DataSourceValue]

    @property
    def value(self) -> Any | None:
        res = self._signal()

        if res is None or not isinstance(res, dict):
            res = self._state.get("default")
            if res is None:
                return None
        
        res_type = res.get("type")
        if res_type == "ldata":
            path = res.get('path')
            if (
                path is None
                or not isinstance(path, str)
                or not path.startswith("latch://")
            ):
                return None
            lpath = LPath(path)

            if path.endswith(".csv"):
                return pd.read_csv(lpath.download())
            elif path.endswith(".xlsx"):
                return pd.read_excel(lpath.download())
            elif path.endswith(".tsv"):
                return pd.read_csv(lpath.download(), sep="\t")
            else:
                return None

        if res_type == "dataframe":
            df_id = res.get('df_id')
            if df_id is None:
                return None

            g = _inject.kernel.k_globals.get_signal(df_id)
            if g is None:
                return None

            return g()

        if res_type == "registry":
            table_id = res.get('table_id')
            if table_id is None:
                return None

            return Table(id=table_id).get_dataframe()


def w_tabular_datasource_picker(
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
