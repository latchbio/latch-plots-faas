from typing import Any

from .. import _inject
from ..reactive import Signal, ctx


def use_state_key(key: str | None = None) -> str:
    if key is not None:
        if "/" in key:
            raise ValueError("'/' is not allowed in widget keys")

        return key

    assert ctx.cur_comp is not None

    res = ctx.cur_comp.widget_state_idx
    ctx.cur_comp.widget_state_idx += 1

    return str(res)


def use_value_signal(key: str) -> Signal[Any]:
    return _inject.kernel.get_widget_value(key)
