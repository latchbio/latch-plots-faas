from typing import TypeAlias

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
