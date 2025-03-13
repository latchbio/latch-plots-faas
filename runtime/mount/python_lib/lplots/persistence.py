from typing import TypedDict


class SerializedSignal(TypedDict):
    value: bytes
    name: str
    listeners: list[str]
    error_msg: None | str
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
