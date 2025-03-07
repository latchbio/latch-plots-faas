from dataclasses import dataclass, field
from datetime import UTC, datetime
from typing import Any, Literal, TypedDict

from ..reactive import Signal
from ..utils.nothing import Nothing
from . import _emit, _state


class ButtonWidgetSignalValue(TypedDict):
    clicked: str
    last_clicked: str


class ButtonWidgetState(_emit.WidgetState[Literal["button"], str]):
    label: str
    readonly: bool
    default: ButtonWidgetSignalValue


datetime_format_string = "%Y-%m-%dT%H:%M:%S.%fZ"


def _iso_string_to_datetime(iso_string: str) -> datetime:
    return datetime.strptime(iso_string, datetime_format_string).replace(tzinfo=UTC)


def parse_iso_strings(data: Any) -> tuple[datetime, datetime] | None:
    if not isinstance(data, dict):
        return None

    required_keys = ButtonWidgetSignalValue.__annotations__.keys()

    if not all(key in data for key in required_keys):
        return None

    try:
        clicked = _iso_string_to_datetime(data["clicked"])
        last_clicked = _iso_string_to_datetime(data["last_clicked"])
        return (clicked, last_clicked)
    except ValueError:
        return None


@dataclass(kw_only=True)
class ButtonWidget:
    _key: str
    _state: ButtonWidgetState
    _signal: Signal[ButtonWidgetSignalValue]

    _last_clicked_ref: None | datetime = field(default=None, repr=False)

    def _value(self, val: ButtonWidgetSignalValue | Nothing) -> bool:
        if val is Nothing.x:
            return False

        parsed = parse_iso_strings(val)
        if parsed is None:
            return False

        clicked, last_clicked = parsed

        if self._last_clicked_ref is None:
            self._last_clicked_ref = last_clicked

        if clicked > self._last_clicked_ref:
            val["last_clicked"] = str(clicked)
            return True

        return False

    @property
    def value(self) -> bool:
        val = self._signal()
        return self._value(val)

    @property
    def sample(self) -> bool:
        res = self._signal.sample()
        return self._value(res)


def w_button(
    *,
    key: str | None = None,
    label: str,
    default: None | ButtonWidgetSignalValue = None,
    readonly: bool = False,
) -> ButtonWidget:
    key = _state.use_state_key(key=key)

    if default is None:
        now = str(datetime.now(UTC).strftime(datetime_format_string))
        default = {"clicked": now, "last_clicked": now}

    res = ButtonWidget(
        _key=key,
        _state={
            "type": "button",
            "label": label,
            "default": default,
            "readonly": readonly,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
