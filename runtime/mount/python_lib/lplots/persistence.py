import subprocess
import sys
from base64 import b64decode, b64encode
from typing import Generic, TypedDict

# todo(kenny): remove when we add uv
try:
    from dill import PicklingError, UnpicklingError, dumps, loads
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "dill"])
    from dill import PicklingError, UnpicklingError, dumps, loads

from .widgets import _emit


class SerializedSignal(TypedDict):
    value: str
    name: str
    listeners: list[str]
    dump_error_msg: str | None
    load_error_msg: str | None
    id: str


class SerializedWidget(TypedDict, Generic[_emit.WidgetType, _emit.WidgetValue]):
    signal_id: str
    state: _emit.WidgetState[_emit.WidgetType, _emit.WidgetValue]
    key: str
    _is_plots_faas_widget: bool


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


class UnUnserialSymbol:
    ...


unable_to_unserialize_symbol: UnUnserialSymbol = UnUnserialSymbol()

MAX_SHORT_VAL_LEN = MAX_REPR_LEN = 100


def safe_serialize_obj(val: object, short: bool = False) -> (str, str | None):
    try:
        s_val = dumps(val)
        error_msg = None
    except (PicklingError, TypeError) as e:
        s_val = dumps(unserial_symbol)
        error_msg = f"Failed to pickle: {e}"
    if short:
        s_val = s_val[:MAX_SHORT_VAL_LEN]
    s_val = b64encode(s_val).decode("utf-8")
    return s_val, error_msg


def safe_unserialize_obj(
    s_val: str,
) -> (object | UnUnserialSymbol, str | None):
    try:
        val = loads(b64decode(s_val.encode("utf-8")))
        error_msg = None
    except UnpicklingError as e:
        error_msg = f"Failed to unpickle: {e}"
        val = unable_to_unserialize_symbol
    return val, error_msg


def small_repr(val: object) -> str:
    val_repr = repr(val)
    if len(val_repr) <= MAX_REPR_LEN:
        return val_repr
    return val_repr[:MAX_REPR_LEN]
