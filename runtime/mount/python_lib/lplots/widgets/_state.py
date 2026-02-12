from .. import _inject
from ..reactive import Signal, get_rctx


def use_state_key(key: str | None = None) -> str:
    if key is not None:
        if "/" in key:
            raise ValueError("'/' is not allowed in widget keys")

        return key

    ctx = get_rctx()
    assert ctx.cur_comp is not None

    res = ctx.cur_comp.widget_state_idx
    ctx.cur_comp.widget_state_idx += 1

    return str(res)


def use_value_signal(key: str) -> Signal[object]:
    return _inject.kernel.get_widget_value(key)


def submit_widget_state() -> None:
    _inject.kernel.submit_widget_state()
