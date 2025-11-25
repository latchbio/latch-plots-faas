from dataclasses import dataclass, field
from datetime import UTC, datetime
from typing import Literal, NotRequired, TypedDict

from ..persistence import SerializedWidget
from ..reactive import Signal
from . import _emit, _state, widget


class ButtonAppearance(TypedDict, total=False):
    variant: Literal["primary", "secondary", "text", "link"]


class ButtonWidgetSignalValue(TypedDict):
    clicked: str
    last_clicked: str


button_type: Literal["button"] = "button"


class ButtonWidgetState(_emit.WidgetState[button_type, str]):
    label: str
    readonly: bool
    default: ButtonWidgetSignalValue
    appearance: NotRequired[ButtonAppearance | None]


class SerializedButtonWidget(SerializedWidget):
    state: ButtonWidgetState  # type: ignore[override]


def parse_iso_strings(data: object) -> tuple[datetime, datetime] | None:
    if not isinstance(data, dict):
        return None

    required_keys = ButtonWidgetSignalValue.__required_keys__

    if not all(key in data for key in required_keys):
        return None

    try:
        clicked = datetime.fromisoformat(data["clicked"])
        last_clicked = datetime.fromisoformat(data["last_clicked"])
        return (clicked, last_clicked)
    except ValueError:
        return None


@dataclass(kw_only=True, frozen=True)
class ButtonWidget(widget.BaseWidget):
    _key: str
    _state: ButtonWidgetState
    _signal: Signal[object | ButtonWidgetSignalValue]

    _last_clicked_ref: list[None | datetime] = field(
        default_factory=lambda: [None], repr=False
    )

    @property
    def value(self) -> bool:
        res = self._signal()

        if not self._signal._ui_update:
            return False

        if not isinstance(res, dict) or not all(
            key in res for key in ButtonWidgetSignalValue.__annotations__
        ):
            return False

        parsed = parse_iso_strings(res)
        if parsed is None:
            return False
        clicked, last_clicked = parsed
        if self._last_clicked_ref[0] is None:
            self._last_clicked_ref[0] = last_clicked

        return clicked > self._last_clicked_ref[0]

    def serialize(self) -> "SerializedButtonWidget":  # type: ignore[override]
        return SerializedButtonWidget(
            signal_id=self._signal.id,
            state=self._state,
            key=self._key,
            _is_plots_faas_widget=True,
        )

    @classmethod
    def load(  # type: ignore[override]
        cls, s_widget: SerializedButtonWidget, widget_sigs: dict[str, Signal]
    ) -> "ButtonWidget":
        sig = widget_sigs[s_widget["signal_id"]]

        return cls(_signal=sig, _state=s_widget["state"], _key=s_widget["key"])


_emit.widget_registry[button_type] = ButtonWidget


def w_button(
    *,
    key: str | None = None,
    label: str,
    default: ButtonWidgetSignalValue | None = None,
    readonly: bool = False,
    appearance: ButtonAppearance | None = None,
) -> ButtonWidget:
    key = _state.use_state_key(key=key)

    if default is None:
        now = datetime.now(UTC).isoformat()
        default = {"clicked": now, "last_clicked": now}

    res = ButtonWidget(
        _key=key,
        _state={
            "type": button_type,
            "label": label,
            "default": default,
            "readonly": readonly,
            "appearance": appearance,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
