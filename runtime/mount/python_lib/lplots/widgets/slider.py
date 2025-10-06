from dataclasses import dataclass
from typing import Literal, NotRequired, TypedDict

from ..reactive import Signal
from . import _emit, _state, widget
from .shared import FormInputAppearance

tool_tip_formatter_type = Literal[
    "number", "percentage", "currency", "decimal", "integer", "scientific", "bytes"
]


class _SliderInputBase(TypedDict):
    label: str
    readonly: bool
    marks: NotRequired[dict[str, int | float]] | None
    step: NotRequired[int | float] | None
    min: NotRequired[int | float] | None
    max: NotRequired[int | float] | None
    tooltip_formatter: NotRequired[tool_tip_formatter_type] | None
    appearance: NotRequired[FormInputAppearance | None]


number_slider_widget_type: Literal["number_slider_input"] = "number_slider_input"


class NumberSliderInputWidgetState(
    _SliderInputBase, _emit.WidgetState[number_slider_widget_type, str]
):
    default: NotRequired[int | float | None]


@dataclass(frozen=True, kw_only=True)
class NumberSliderInputWidget(widget.BaseWidget):
    _key: str
    _state: NumberSliderInputWidgetState
    _signal: Signal[int | float]

    def _value(self, val: int | float | None) -> int | float:
        if val is None:
            default = self._state.get("default")
            if default is None:
                return self._state.get("min", 0)

            return default

        return val

    @property
    def value(self) -> int | float:
        res = self._signal()
        return self._value(res)

    def sample(self) -> int | float:
        res = self._signal.sample()
        return self._value(res)


_emit.widget_registry[number_slider_widget_type] = NumberSliderInputWidget


def w_number_slider_input(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    default: int | float | None = None,
    min: int | float | None = None,
    max: int | float | None = None,
    step: int | float | None = None,
    tooltip_formatter: tool_tip_formatter_type | None = None,
    marks: dict[str, int | float] | None = None,
    appearance: FormInputAppearance | None = None,
) -> NumberSliderInputWidget:
    key = _state.use_state_key(key=key)

    res = NumberSliderInputWidget(
        _key=key,
        _state={
            "type": number_slider_widget_type,
            "label": label,
            "default": default,
            "readonly": readonly,
            "min": min,
            "max": max,
            "step": step,
            "tooltip_formatter": tooltip_formatter,
            "marks": marks,
            "appearance": appearance,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res


range_slider_widget_type: Literal["range_slider_input"] = "range_slider_input"

RangeSliderValue = tuple[int, int] | tuple[float, float]


class RangeSliderInputWidgetState(
    _SliderInputBase, _emit.WidgetState[range_slider_widget_type, str]
):
    default: NotRequired[RangeSliderValue | None]


@dataclass(frozen=True, kw_only=True)
class RangeSliderInputWidget(widget.BaseWidget):
    _key: str
    _state: RangeSliderInputWidgetState
    _signal: Signal[RangeSliderValue]

    def _value(self, val: RangeSliderValue | None) -> RangeSliderValue:
        if val is None:
            default = self._state.get("default")
            if default is None:
                min_value = self._state.get("min", 0)
                max_value = self._state.get("max", 100)
                return (min_value, max_value)

            return default

        return val

    @property
    def value(self) -> RangeSliderValue:
        res = self._signal()
        return self._value(res)

    def sample(self) -> RangeSliderValue:
        res = self._signal.sample()
        return self._value(res)


_emit.widget_registry[range_slider_widget_type] = RangeSliderInputWidget


def w_range_slider_input(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    default: RangeSliderValue | None = None,
    min: int | float | None = None,
    max: int | float | None = None,
    step: int | float | None = None,
    tooltip_formatter: tool_tip_formatter_type | None = None,
    marks: dict[str, int | float] | None = None,
    appearance: FormInputAppearance | None = None,
) -> RangeSliderInputWidget:
    key = _state.use_state_key(key=key)

    res = RangeSliderInputWidget(
        _key=key,
        _state={
            "type": range_slider_widget_type,
            "label": label,
            "default": default,
            "readonly": readonly,
            "min": min,
            "max": max,
            "step": step,
            "tooltip_formatter": tooltip_formatter,
            "marks": marks,
            "appearance": appearance,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res
