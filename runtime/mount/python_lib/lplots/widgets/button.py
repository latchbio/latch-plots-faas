import asyncio
from dataclasses import dataclass
from datetime import UTC, datetime
from typing import Any, Literal, TypedDict

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
    _trigger_signal: Signal[object | int]
    _last_observed_generation: int = 0

    @property
    def value(self) -> bool:
        trg_sig_gen = self._trigger_signal()
        if (
            not isinstance(trg_sig_gen, int)
            or self._last_observed_generation > trg_sig_gen
        ):
            return False

        self._last_observed_generation = trg_sig_gen
        return True

    def _update(self) -> None:
        res = self._signal()

        if not isinstance(res, dict) or not all(
            key in res for key in ButtonWidgetSignalValue.__annotations__
        ):
            return

        parsed = parse_iso_strings(res)
        if parsed is None:
            return

        clicked, last_clicked = parsed

        if clicked <= last_clicked:
            return

        self._signal({**res, "last_clicked": str(clicked)})

        trg_sig_gen = self._trigger_signal.sample()
        if not isinstance(trg_sig_gen, int):
            trg_sig_gen = 0

        self._trigger_signal(trg_sig_gen + 1)

    async def _create_update_node(self) -> None:
        await ctx.run(self._update)


def w_button(
    *,
    key: str | None = None,
    label: str,
    default: None | ButtonWidgetSignalValue = None,
    readonly: bool = False,
) -> ButtonWidget:
    key = _state.use_state_key(key=key)
    trigger_key = _state.use_state_key(key=f"{key}#trigger")

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
        _trigger_signal=_state.use_value_signal(key=trigger_key),
    )
    _emit.emit_widget(key, res._state)

    # todo(rteqs): this can deadlock. either we make w_button async or figure something out
    asyncio.run_coroutine_threadsafe(
        res._create_update_node(), asyncio.get_running_loop()
    )

    return res
