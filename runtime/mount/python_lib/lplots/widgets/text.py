from dataclasses import dataclass
from typing import Literal, NotRequired, TypedDict

from ..reactive import Signal
from ..utils.nothing import Nothing
from . import _emit, _state
from .shared import FormInputAppearance


class TextInputWidgetState(_emit.WidgetState[Literal["text_input"], str]):
    label: str
    readonly: bool
    default: NotRequired[str | None]
    appearance: NotRequired[FormInputAppearance | None]


@dataclass(frozen=True, kw_only=True)
class TextInputWidget:
    _key: str
    _state: TextInputWidgetState
    _signal: Signal[str]

    @property
    def value(self) -> str:
        res = self._signal()
        if res is Nothing.x or not isinstance(res, str):
            res = self._state.get("default")

            if res is None:
                res = ""

        return res


def w_text_input(
    *,
    key: str | None = None,
    label: str,
    readonly: bool = False,
    default: str | None = None,
    appearance: FormInputAppearance | None = None,
) -> TextInputWidget:
    key = _state.use_state_key(key=key)

    res = TextInputWidget(
        _key=key,
        _state={
            "type": "text_input",
            "label": label,
            "default": default,
            "readonly": readonly,
            "appearance": appearance,
        },
        _signal=_state.use_value_signal(key=key),
    )
    _emit.emit_widget(key, res._state)

    return res


class TextOutputAppearance(TypedDict, total=False):
    message_box: Literal["danger", "info", "success", "warning", "primary", "neutral"]


class TextOutputWidgetState(_emit.WidgetState[Literal["text_output"], None]):
    content: str
    appearance: NotRequired[TextOutputAppearance | None]


def w_text_output(
    *,
    key: str | None = None,
    content: str,
    appearance: TextOutputAppearance | None = None,
) -> None:
    key = _state.use_state_key(key=key)

    w: TextOutputWidgetState = {
        "type": "text_output",
        "content": content,
        "appearance": appearance,
    }
    _emit.emit_widget(key, w)
