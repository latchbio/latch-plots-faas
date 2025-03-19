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

    _clicked_ref: None | datetime = field(default=None, repr=False)
    _last_clicked_ref: None | datetime = field(default=None, repr=False)

    @property
    def value(self) -> bool:
        self._lambda_signal()

        if self._clicked_ref is None or self._last_clicked_ref is None:
            return False

        return self._clicked_ref > self._last_clicked_ref

    def h(self) -> None:
        print("DEBUG: _helper")
        self._lambda_signal(None)

        res = self._signal()

        if not isinstance(res, dict) or not all(
            key in res for key in ButtonWidgetSignalValue.__annotations__
        ):
            return

        parsed = parse_iso_strings(res)
        if parsed is None:
            return

        clicked, last_clicked = parsed

        if self._last_clicked_ref is None:
            self._last_clicked_ref = last_clicked

        self._clicked_ref = clicked

        if clicked > self._last_clicked_ref:
            res["last_clicked"] = str(clicked)

        print("DEBUG: _helper done")
        return

    async def _helper(self) -> None:
        async with ctx.transaction:
            await ctx.run(self.h)


def w_button(
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

    lambda_signal = _state.use_value_signal(key=lambda_key)
    res = ButtonWidget(
        _key=key,
        _state={
            "type": "button",
            "label": label,
            "default": default,
            "readonly": readonly,
        },
        _signal=_state.use_value_signal(key=key),
        _lambda_signal=lambda_signal,
    )
    _emit.emit_widget(key, res._state)

    if lambda_signal.sample() is Nothing.x:
        asyncio.run_coroutine_threadsafe(
            res._helper(), asyncio.get_running_loop()
        ).done()

    return res
