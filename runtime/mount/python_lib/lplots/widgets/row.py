from typing import Literal

from . import _emit, _state


class RowWidgetState(_emit.WidgetState[Literal["row"], None]):
    items: list[str]


def w_row(
    *,
    key: str | None = None,
    items: set[object],
) -> None:
    key = _state.use_state_key(key=key)

    # todo(manske): have better type checking here for widgets
    assert all(hasattr(i, "_key") for i in items)
    item_keys = [i._key for i in items]

    w: RowWidgetState = {
        "type": "row",
        "items": set(item_keys),
    }

    _emit.emit_widget(key, w)
