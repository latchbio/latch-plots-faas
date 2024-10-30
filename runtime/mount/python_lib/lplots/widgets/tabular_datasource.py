from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from typing import Any, Literal, NotRequired, Union

import pandas as pd
from latch.ldata.path import LPath
from latch.registry.table import Table

from runtime.mount.lplots.reactive import Signal

from .. import _inject
from . import _emit, _state
from .select import Select, w_select
from .shared import FormInputAppearance

DataSourceType = Literal["ldata", "dataframe", "registry"]


@dataclass
class LDataDataSource:
    path: str
    type: Literal["ldata"]


@dataclass
class DataFrameDataSource:
    df_id: str
    type: Literal["dataframe"]


@dataclass
class RegistryDataSource:
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
        if res is None:
            return None
        elif res.type == "ldata":
            path = res.path
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
        elif res.type == "dataframe":
            df_id = res.df_id
            if df_id is None:
                return None

            g = _inject.kernel.k_globals.get_signal(df_id)
            if g is None:
                return None

            return g()
        elif res.type == "registry":
            table_id = res.table_id
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
