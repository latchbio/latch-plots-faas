import asyncio
from dataclasses import dataclass, field
from datetime import UTC, datetime
from typing import Any, Literal, TypedDict

from lplots.utils.nothing import Nothing

from ..reactive import Signal, ctx
from . import _emit, _state


class ButtonWidgetSignalValue(TypedDict):
    clicked: str
    last_clicked: str


class ButtonWidgetState(_emit.WidgetState[Literal["button"], str]):
    label: str
    readonly: bool
    default: ButtonWidgetSignalValue


def parse_iso_strings(data: Any) -> tuple[datetime, datetime] | None:
    if not isinstance(data, dict):
        return None

    required_keys = ButtonWidgetSignalValue.__annotations__.keys()

    if not all(key in data for key in required_keys):
        return None

    try:
        clicked = datetime.fromisoformat(data["clicked"])
        last_clicked = datetime.fromisoformat(data["last_clicked"])
        return (clicked, last_clicked)
    except ValueError:
        return None


@dataclass(kw_only=True)
class ButtonWidget:
    _key: str
    _state: ButtonWidgetState
    _signal: Signal[object | ButtonWidgetSignalValue]
    _lambda_signal: Signal[object]

    _last_clicked_ref: None | datetime = field(default=None, repr=False)

    @property
    def value(self) -> bool:
        self._lambda_signal()
        res = self._signal.sample()

        if not isinstance(res, dict) or not all(
            key in res for key in ButtonWidgetSignalValue.__annotations__
        ):
            return False

        parsed = parse_iso_strings(res)
        if parsed is None:
            return False

        clicked, last_clicked = parsed

        if self._last_clicked_ref is None:
            self._last_clicked_ref = last_clicked

        if clicked > self._last_clicked_ref:
            self._signal({"clicked": str(clicked), "last_clicked": str(last_clicked)})
            return True

        return False

    def _helper(self) -> None:
        self._signal()
        self._lambda_signal(Nothing.x)


async def w_button(
    *,
    key: str | None = None,
    label: str,
    default: None | ButtonWidgetSignalValue = None,
    readonly: bool = False,
) -> ButtonWidget:
    key = _state.use_state_key(key=key)
    lambda_key = _state.use_state_key(key=f"lambda_{key}")

    if default is None:
        now = datetime.now(UTC).isoformat()
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
        _lambda_signal=_state.use_value_signal(key=lambda_key),
    )
    _emit.emit_widget(key, res._state)

    await ctx.run(res._helper)
    return res
