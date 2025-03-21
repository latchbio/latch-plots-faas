from base64 import b64decode, b64encode
from typing import TypedDict

from dill import dumps, loads


class SerializedSignal(TypedDict):
    value: str
    name: str
    listeners: list[str]
    error_msg: str | None
    id: str


class SerializedNode(TypedDict):
    code: str
    stale: bool
    signals: list[str]
    cell_id: str | None
    name: str | None
    parent: str | None
    id: str
    widget_state_idx: int
    widget_states: list[dict]


unserial_symbol = "<<UNSERIALIZABLE>>"

MAX_SHORT_VAL_LEN = MAX_REPR_LEN = 1000


def safe_serialize_obj(val: object, short: bool = False) -> (str, str | None):
    try:
        s_val = dumps(val)
        error_msg = None
    except Exception as e:
        s_val = dumps(unserial_symbol)
        error_msg = f"Failed to pickle: {e}"
    if short:
        s_val = s_val[:MAX_SHORT_VAL_LEN]
    s_val = b64encode(s_val).decode("utf-8")
    return s_val, error_msg


def safe_unserialize_obj(s_val: str) -> object | None:
    try:
        val = loads(b64decode(s_val.encode("utf-8")))
    except Exception as e:
        # todo(kenny): is it OK to collapse with actual None
        val = None
    return val


def small_repr(val: object) -> str:
    val_repr = repr(val)
    if len(val_repr) <= MAX_REPR_LEN:
        return val_repr
    return val_repr[:MAX_REPR_LEN]
