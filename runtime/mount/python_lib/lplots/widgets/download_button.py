from dataclasses import dataclass
from typing import Literal, NotRequired

from . import _emit, _state, widget
from .button import ButtonAppearance

download_button_type: Literal["download_button"] = "download_button"


class DownloadButtonWidgetState(_emit.WidgetState[download_button_type, str]):
    label: str
    source: str
    appearance: NotRequired[ButtonAppearance | None]


@dataclass(frozen=True, kw_only=True)
class DownloadButtonWidget(widget.BaseWidget):
    _key: str
    _state: DownloadButtonWidgetState


_emit.widget_registry[download_button_type] = DownloadButtonWidget


def w_download_button(
    *,
    key: str | None = None,
    label: str,
    source: str,
    appearance: ButtonAppearance | None = None,
) -> DownloadButtonWidget:
    key = _state.use_state_key(key=key)

    res = DownloadButtonWidget(
        _key=key,
        _state={
            "type": download_button_type,
            "label": label,
            "source": source,
            "appearance": appearance,
        },
    )
    _emit.emit_widget(key, res._state)

    return res
