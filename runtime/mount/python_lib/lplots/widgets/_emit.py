from typing import TYPE_CHECKING, Generic, NotRequired, TypedDict, TypeVar

from .. import _inject

if TYPE_CHECKING:
    from ..reactive import Signal

WidgetType = TypeVar("WidgetType", bound=str)
WidgetValue = TypeVar("WidgetValue")


class WidgetState(TypedDict, Generic[WidgetType, WidgetValue]):
    type: WidgetType
    value: NotRequired["Signal[WidgetValue]"]


def emit_widget(key: str, data: WidgetState):
    _inject.kernel.emit_widget(key, data)
