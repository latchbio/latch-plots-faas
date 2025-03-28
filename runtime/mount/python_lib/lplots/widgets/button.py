import asyncio
from dataclasses import dataclass
from datetime import UTC, datetime
from typing import Literal, TypedDict

from ..persistence import SerializedWidget
from ..reactive import Signal, ctx
from . import _emit, _state, widget


class ButtonWidgetSignalValue(TypedDict):
    clicked: str
    last_clicked: str


class ButtonWidgetState(_emit.WidgetState[Literal["button"], str]):
    label: str
    readonly: bool
    default: ButtonWidgetSignalValue


class SerializedButtonWidget(SerializedWidget):
    trigger_signal_id: str
    # state: _emit.WidgetState[Literal["button"], str]
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


@dataclass(frozen=True, kw_only=True)
class ButtonWidget(widget.BaseWidget):
    _key: str
    _state: ButtonWidgetState
    _signal: Signal[object | ButtonWidgetSignalValue]
    _trigger_signal: Signal[object]

    @property
    def value(self) -> bool:
        self._trigger_signal()
        return self._trigger_signal.id in ctx.prev_updated_signals

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
        self._trigger_signal(None)

    async def _create_update_node(self) -> None:
        await ctx.run(self._update)

    def serialize(self) -> "SerializedButtonWidget":  # type: ignore[override]
        return SerializedButtonWidget(
            signal_id=self._signal.id,
            trigger_signal_id=self._trigger_signal.id,
            state=self._state,
            key=self._key,
        )

    @classmethod
    def load(  # type: ignore[override]
        cls, s_widget: SerializedButtonWidget, widget_sigs: dict[str, Signal]
    ) -> "ButtonWidget":
        sig = widget_sigs[s_widget["signal_id"]]
        trigger_sig = widget_sigs[s_widget["trigger_signal_id"]]

        return cls(
            _signal=sig,
            _trigger_signal=trigger_sig,
            _state=s_widget["state"],
            _key=s_widget["key"],
        )


def w_button(
    *,
    key: str | None = None,
    label: str,
    default: ButtonWidgetSignalValue | None = None,
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
