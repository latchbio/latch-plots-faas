from abc import ABC
from dataclasses import dataclass, field

from ..persistence import SerializedWidget
from ..reactive import Signal
from . import _emit


@dataclass(frozen=True, kw_only=True)
class BaseWidget(ABC):
    _key: str
    _state: _emit.WidgetState

    # todo(maximsmol): fix typing here and make signal fields optional without breaking serialization somehow
    _signal: Signal = field(default_factory=lambda: Signal(None))
    _has_signal: bool = True

    def serialize(self) -> SerializedWidget:
        return SerializedWidget(
            signal_id=self._signal.id,
            state=self._state,
            key=self._key,
            _is_plots_faas_widget=True,
        )

    @classmethod
    def load(
        cls, s_widget: SerializedWidget, widget_sigs: dict[str, Signal]
    ) -> "BaseWidget":
        # todo(kenny): unsure if we want to throw here if not in dict
        sig = widget_sigs[s_widget["signal_id"]]
        return cls(_signal=sig, _state=s_widget["state"], _key=s_widget["key"])


def load_widget_helper(
    s_widget: SerializedWidget, widget_sigs: dict[str, Signal]
) -> "BaseWidget":
    # todo(kenny): unsure if we want to throw here if not in dict
    w_cls = _emit.widget_registry.get(s_widget["state"]["type"])
    assert w_cls is not None, f"{s_widget['state']['type']} registry is None"
    return w_cls.load(s_widget, widget_sigs)
