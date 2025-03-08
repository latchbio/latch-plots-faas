from dataclasses import dataclass
from typing import Any

from .. import _inject
from . import _emit, _state
from .select import Select, w_select
from .shared import FormInputAppearance


@dataclass(frozen=True, kw_only=True)
class DataframePicker:
    _select: Select

    @property
    def key(self) -> str | None:
        res = self._select.value
        assert isinstance(res, str)
        return res

    def _value(self, key: str | None) -> Any | None:
        if key is None:
            return None

        g = _inject.kernel.k_globals.get_signal(key)
        if g is None:
            return None

        return g()

    @property
    def value(self) -> Any | None:
        key = self._select.value
        assert isinstance(key, str)
        return self._value(key)

    def sample(self) -> Any | None:
        key = self._select.sample
        assert isinstance(key, str)
        return self._value(key)


def w_dataframe_picker(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    appearance: FormInputAppearance | None = None,
    required: bool = False,
) -> DataframePicker:
    key = _state.use_state_key(key=key)

    dfs = _inject.kernel.k_globals.dataframes()
    res = w_select(
        key=key,
        label=label,
        readonly=readonly,
        appearance=appearance,
        options=sorted(dfs),
        required=required,
    )
    _emit.emit_widget(key, res._state)

    return DataframePicker(_select=res)
