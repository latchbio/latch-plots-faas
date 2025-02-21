from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Literal, Tuple, TypedDict

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


def _iso_string_to_datetime(iso_string: str) -> datetime:
    return datetime.strptime(iso_string, "%Y-%m-%dT%H:%M:%S.%fZ").replace(
        tzinfo=timezone.utc
    )


def parse_iso_strings(data: Any) -> None | Tuple[datetime, datetime]:
    if not isinstance(data, dict):
        return

    required_keys = ButtonWidgetSignalValue.__annotations__.keys()

    if not all(key in data for key in required_keys):
        return

    try:
        clicked = _iso_string_to_datetime(data["clicked"])
        last_clicked = _iso_string_to_datetime(data["last_clicked"])
        return (clicked, last_clicked)
    except ValueError:
        return


@dataclass(kw_only=True)
class ButtonWidget:
    _key: str
    _state: ButtonWidgetState
    _signal: Signal[ButtonWidgetSignalValue]

    _last_clicked_ref: None | datetime = field(default=None, repr=False)

    @property
    def value(self) -> bool:
        res = self._signal()

        if res is Nothing.x:
            return False

        parsed = parse_iso_strings(res)
        if parsed is None:
            return False
        clicked, last_clicked = parsed

        if self._last_clicked_ref is None:
            self._last_clicked_ref = last_clicked

        if clicked > self._last_clicked_ref:
            res["last_clicked"] = str(clicked)
            return True

        return False


def w_button(
    *,
    key: str | None = None,
    label: str,
    default: None | ButtonWidgetSignalValue = None,
    readonly: bool = False,
) -> ButtonWidget:
    key = _state.use_state_key(key=key)

    if default is None:
        now = str(datetime.now(timezone.utc))
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
