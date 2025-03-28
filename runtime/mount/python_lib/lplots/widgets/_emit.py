from typing import TYPE_CHECKING, Generic, NotRequired, TypedDict, TypeVar

from .. import _inject
from . import widget

if TYPE_CHECKING:
    from ..reactive import Signal

WidgetType = TypeVar("WidgetType", bound=str)
WidgetValue = TypeVar("WidgetValue")

widget_registry: dict[str, widget.BaseWidget] = {}


class WidgetState(TypedDict, Generic[WidgetType, WidgetValue]):
    type: WidgetType
    value: NotRequired["Signal[WidgetValue]"]


def emit_widget(key: str, data: WidgetState[WidgetType, WidgetValue]) -> None:
    _inject.kernel.emit_widget(key, data)
