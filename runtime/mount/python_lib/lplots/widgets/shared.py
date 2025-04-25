from typing import TypeAlias, TypedDict

from .ann_data import AnnData
from .checkbox import CheckboxWidget
from .dataframe import DataframePicker
from .datasource import TabularDatasourcePicker
from .ldata import LDataPicker
from .multiselect import MultiSelect
from .plot import Plot
from .radio import RadioGroups
from .registry import RegistryTablePicker
from .select import Select
from .table import Table
from .text import TextInputWidget


class FormInputAppearance(TypedDict, total=False):
    detail: str | None
    placeholder: str | None
    help_text: str | None
    error_text: str | None
    description: str | None


InteractiveWidget: TypeAlias = (
    CheckboxWidget
    | LDataPicker
    | Select
    | RadioGroups
    | MultiSelect
    | DataframePicker
    | TextInputWidget
    | RegistryTablePicker
    | TabularDatasourcePicker
    | Plot
    | Table
    | AnnData
)
