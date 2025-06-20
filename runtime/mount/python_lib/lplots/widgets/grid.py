import types
from dataclasses import dataclass
from typing import Literal, Self

from . import _emit, _state, widget
from .widget import BaseWidget

grid_widget_type: Literal["grid"] = "grid"


# todo(manske): allow for custom css to define layout
@dataclass(frozen=True, kw_only=True)
class GridItem:
    item: str
    col_span: int
    row_span: int


class GridWidgetState(_emit.WidgetState[grid_widget_type, None]):
    grid_items: list[GridItem]
    template: str


@dataclass(frozen=True, kw_only=True)
class Grid(widget.BaseWidget):
    _key: str
    _state: GridWidgetState

    _has_signal = False

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: types.TracebackType | None,
    ) -> None:
        pass

    def add(
        self,
        item: BaseWidget,
        col_span: int = 1,
        row_span: int = 1,
    ) -> None:
        grid_item = GridItem(
            item=item._key,
            col_span=col_span,
            row_span=row_span,
        )

        return self._state["grid_items"].append(grid_item)


# todo(manske): allow for custom css to define layout
def w_grid(
    *,
    key: str | None = None,
    columns: int,
    rows: int | None = None,
) -> Grid:
    key = _state.use_state_key(key=key)

    column_template = f"repeat({columns}, 1fr)"

    row_template = "none"
    if rows is not None:
        row_template = f"repeat({rows}, auto)"

    template = f"{row_template} / {column_template}"

    res = Grid(
        _key=key,
        _state={
            "type": "grid",
            "grid_items": [],
            "template": template,
        },
    )

    _emit.emit_widget(key, res._state)
    return res
