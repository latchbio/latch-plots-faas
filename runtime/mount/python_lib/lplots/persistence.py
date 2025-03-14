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


unserial_symbol = "<<UNSERIALIZABLE>>"


def safe_serialize_obj(val: object) -> (str, str | None):
    try:
        s_val = dumps(val)
        error_msg = None
    except Exception as e:
        s_val = dumps(unserial_symbol)
        error_msg = f"Failed to pickle: {e}"
    s_val = b64encode(s_val).decode("utf-8")
    return s_val, error_msg


def safe_unserialize_obj(s_val: str) -> object | None:
    try:
        val = loads(b64decode(s_val.encode("utf-8")))
    except Exception as e:
        # todo(kenny): is it OK to collapse with actual None
        val = None
    return val
