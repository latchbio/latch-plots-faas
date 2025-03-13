from typing import TypedDict


class SerializedSignal(TypedDict):
    value: bytes
    name: str
    listeners: list[str]
    error_msg: None | str
    id: str


class SerializedNode(TypedDict):
    code: str
    stale: str
    signals: list[str]
    cell_id: str
    name: str
    parent: str
    id: str


unserial_symbol = "<<UNSERIALIZABLE>>"
