from abc import ABC
from dataclasses import dataclass

from ..persistence import SerializedWidget
from ..reactive import Signal
from . import _emit


@dataclass(frozen=True, kw_only=True)
class BaseWidget(ABC):
    _key: str
    _state: _emit.WidgetState
    _signal: Signal

    def serialize(self) -> SerializedWidget:
        return SerializedWidget(
            signal_id=self._signal.id, state=self._state, key=self._key
        )

    @classmethod
    def load(
        cls, s_widget: SerializedWidget, widget_sigs: dict[str, Signal]
    ) -> "BaseWidget":
        # todo(kenny): think we want to throw here if not here but unsure
        sig = widget_sigs[s_widget["signal_id"]]
        return cls(_signal=sig, _state=s_widget["state"], _key=s_widget["key"])
