import inspect
import sys
import uuid
from collections.abc import AsyncGenerator, Awaitable, Callable
from contextlib import asynccontextmanager
from dataclasses import dataclass, field
from traceback import print_exc
from typing import Any, Generic, Self, TextIO, TypeAlias, TypeVar, overload

from . import _inject
from .persistence import (
    SerializedNode,
    SerializedSignal,
    safe_serialize_obj,
    safe_unserialize_obj,
    unable_to_unserialize_symbol,
)
from .utils.nothing import Nothing
from .widgets._emit import WidgetState

T = TypeVar("T")
R = TypeVar("R")
Computation: TypeAlias = Callable[..., R]

live_nodes: dict[str, "Node"] = {}
live_node_ids: set[str] = set()

live_signals: dict[str, "Signal[object]"] = {}
live_signal_ids: set[str] = set()

stub_node_noop = lambda: None


def graphviz() -> None:
    from pathlib import Path

    with Path("graph.dot").open("w", encoding="utf-8") as f:
        f.write("digraph G {\n")

        for n in live_nodes.values():
            n.graphviz(f)

        f.write("}\n")


@dataclass(kw_only=True)
class Node:
    f: Computation[Any]
    code: str
    stale: bool = False
    disposed: bool = False

    parent: Self | None
    children: dict[str, Self] = field(default_factory=dict)

    signals: dict[str, "Signal"] = field(default_factory=dict)

    cell_id: str | None = None
    name: str | None = None
    _id: str | None

    _is_stub: bool = False

    widget_states: dict[str, WidgetState[str, object]] = field(default_factory=dict)
    widget_state_idx: int = 0

    @property
    def id(self) -> str:
        assert self._id is not None
        return self._id

    def __post_init__(self) -> None:
        if self._id is None:
            self._id = str(uuid.uuid4())

        if self.parent is not None:
            self.parent.children[self.id] = self
            self.cell_id = self.parent.cell_id
        elif self.name is None:
            self.name = self.cell_id

        if self.name is None:
            self.name = self.f.__name__

        if self.id in live_node_ids:
            raise ValueError(f"reactive node id is not unique: {self.id!r}")

        live_nodes[self.id] = self
        live_node_ids.add(self.id)

    def serialize(self) -> SerializedNode:
        return SerializedNode(
            code=self.code,
            stale=self.stale,
            signals=list(self.signals.keys()),
            cell_id=self.cell_id,
            name=self.name,
            parent=(self.parent.id if self.parent else None),
            id=self.id,
            widget_state_idx=self.widget_state_idx,
            widget_states=self.widget_states,
        )

    @classmethod
    def load(cls, s_node) -> "Node":
        return cls(
            f=stub_node_noop,
            code=s_node["code"],
            stale=s_node["stale"],
            signals={},
            cell_id=s_node["cell_id"],
            name=s_node["name"],
            parent=None,
            widget_state_idx=s_node["widget_state_idx"],
            widget_states=s_node["widget_states"],
            _id=s_node["id"],
            _is_stub=True,
        )

    def name_path(self) -> str:
        assert self.name is not None
        res = [self.name]

        cur = self.parent
        while cur is not None:
            assert cur.name is not None

            res.append(cur.name)
            cur = cur.parent

        return "/".join(reversed(res))

    def farthest_stale_ancestor(self) -> Self | None:
        res: Self | None = None

        cur = self
        while cur is not None:
            if cur.stale:
                res = cur

            cur = cur.parent

        return res

    def dispose(self) -> None:
        _inject.kernel.on_dispose(self)

        stack = [self]

        while len(stack) > 0:
            cur = stack.pop()
            # print("[@#] disposed", cur)

            assert not cur.disposed
            cur.disposed = True

            if cur.parent is not None:
                del cur.parent.children[cur.id]
            cur.parent = None

            stack.extend(cur.children.values())
            # we let the children delete themselves
            # this is because the topmost node needs to remove itself from a potential parent
            # which is not itself being disposed

            for s in cur.signals.values():
                del s._listeners[cur.id]
            cur.signals = {}

            del live_nodes[cur.id]

    def __repr__(self) -> str:
        stale_mark = "!" if self.stale else ""
        return f"{stale_mark}{self.name}#{self.cell_id}@{self.id}<{self.parent}"

    def debug_state(self, *, no_parent: bool = False) -> dict[str, object]:
        return {
            "repr": repr(self),
            "f": {"qualname": self.f.__qualname__, "name": self.f.__name__},
            "stale": self.stale,
            "disposed": self.disposed,
            "parent": "omitted"
            if no_parent
            else self.parent.debug_state()
            if self.parent is not None
            else None,
            "children": {
                str(k): v.debug_state(no_parent=True) for k, v in self.children.items()
            },
            "signals": {str(k): repr(v) for k, v in self.signals.items()},
            "cell_id": self.cell_id,
            "name": self.name,
            "widget_states": self.widget_states,
            "widget_state_idx": self.widget_state_idx,
        }

    def graphviz(self, f: TextIO) -> None:
        sig_list = [x._name for x in self.signals.values()]

        name = self.f.__name__.replace("<", "&lt;").replace(">", "&gt;")
        label_str = "<BR />".join([f"{name} @ {self.id}", ", ".join(sig_list)])

        f.write(f"{self.id}[label=<{label_str}>];\n")

        if self.parent is not None:
            f.write(f"{self.parent.id} -> {self.id};\n")


@dataclass
class RCtx:
    cur_comp: Node | None = None

    updated_signals: dict[str, "Signal"] = field(default_factory=dict)
    signals_updated_from_code: dict[str, "Signal"] = field(default_factory=dict)
    stale_nodes: dict[str, Node] = field(default_factory=dict)
    prev_updated_signals: dict[str, "Signal"] = field(default_factory=dict)

    in_tx: bool = False

    async def run(
        self,
        f: Callable[..., Awaitable[R]],
        code: str | None = None,
        *,
        _cell_id: str | None = None,
    ) -> R:
        # note(kenny): it is only safe to call without code if the node is a
        # child and can be reconstructed by the parent. Otherwise not possible
        # to serialize when snapshot requested.
        assert code is not None or self.cur_comp is not None

        # note(maximsmol): we want this to happen for non-cell nodes too
        # so it has to be inside `RCtx` which sees every ran node
        # and not just the cell body in the kernel
        if _cell_id is not None:
            await _inject.kernel.set_active_cell(_cell_id)

        async with self.transaction:
            self.cur_comp = Node(
                f=f, parent=self.cur_comp, cell_id=_cell_id, _id=None, code=code
            )

            try:
                if inspect.iscoroutinefunction(f):
                    return await f()
                return f()
            finally:
                self.cur_comp = self.cur_comp.parent

    async def _tick(self) -> None:
        tick_updated_signals = {
            **self.signals_updated_from_code,
            **self.updated_signals,
        }
        self.signals_updated_from_code = {}

        try:
            stack_depth = 1

            f = inspect.currentframe()
            while f is not None:
                f = f.f_back
                stack_depth += 1

            if stack_depth > sys.getrecursionlimit() - 10:
                raise RecursionError("maximum recursion depth exceeded")

            assert self.cur_comp is None
            assert not self.in_tx

            if len(self.updated_signals) == 0:
                return

            for s in self.updated_signals.values():
                s._apply_updates()

            self.prev_updated_signals = self.updated_signals
            self.updated_signals = {}

            to_dispose: dict[str, tuple[Node, Node | None]] = {}
            for n in self.stale_nodes.values():
                if n.disposed:
                    continue

                fsa = n.farthest_stale_ancestor()
                assert fsa is not None

                to_dispose[fsa.id] = (fsa, fsa.parent)

            self.stale_nodes = {}

            for n, _p in to_dispose.values():
                n.dispose()

            if len(to_dispose) > 1:
                async with self.transaction:
                    for n, p in to_dispose.values():
                        self.cur_comp = p

                        try:
                            if n._is_stub:
                                n._is_stub = False
                                assert n.cell_id is not None
                                # reconstruct the function with globals
                                await _inject.kernel.exec(
                                    cell_id=n.cell_id, code=n.code, _from_stub=True
                                )
                            else:
                                await self.run(n.f, n.code, _cell_id=n.cell_id)

                        except Exception:
                            print_exc()
                        finally:
                            self.cur_comp = None

        finally:
            await _inject.kernel.on_tick_finished(tick_updated_signals)
            self.prev_updated_signals = {}
            for sig in live_signals.values():
                sig._ui_update = False

    @property
    @asynccontextmanager
    async def transaction(self) -> AsyncGenerator[None, None]:
        if self.in_tx:
            yield
            return

        try:
            self.in_tx = True
            yield
        finally:
            self.in_tx = False
            await self._tick()


ctx = RCtx()


class Updater(Generic[T]):
    f: Callable[[T], T]


@dataclass(init=False)
class Signal(Generic[T]):
    _value: T
    _name: str
    _id: str

    _updates: list[T | Updater[T]]
    _listeners: dict[str, Node]

    _ui_update: bool = False

    _load_error_msg: str | None

    def __init__(
        self,
        initial: T,
        *,
        name: str | None = None,
        _id: str | None = None,
        _load_error_msg: str | None = None,
    ) -> None:
        if _id is None:
            self._id = str(uuid.uuid4())
        else:
            self._id = _id

        if name is None:
            name = f"Signal@{self.id}"

        self._value = initial
        self._name = name

        self._updates = []
        self._listeners = {}

        if self.id in live_signal_ids:
            raise ValueError(f"signal id is not unique: {self.id!r}")

        live_signals[self.id] = self
        live_signal_ids.add(self.id)

        self._load_error_msg = _load_error_msg

    @property
    def id(self) -> str:
        assert self._id is not None
        return self._id

    def sample(self) -> T:
        return self._value

    def serialize(self, short_val: bool = False) -> SerializedSignal:
        s_val, error_msg = safe_serialize_obj(self._value, short=short_val)
        return SerializedSignal(
            value=s_val,
            name=self._name,
            listeners=list(self._listeners.keys()),
            dump_error_msg="" if error_msg is None else error_msg,
            load_error_msg=self._load_error_msg,
            id=self.id,
        )

    @classmethod
    def load(cls, s_sig) -> "Signal[T]":
        val, error_msg = safe_unserialize_obj(s_sig["value"])

        if val is None or val is unable_to_unserialize_symbol:
            sig = cls(
                Nothing.x,
                name=s_sig["name"],
                _id=s_sig["id"],
                _load_error_msg=error_msg,
            )
        else:
            sig = cls(val, name=s_sig["name"], _id=s_sig["id"])

        return sig

    @overload
    def __call__(self, /) -> T:
        ...

    @overload
    def __call__(self, /, upd: T, *, _ui_update: bool = False) -> None:
        ...

    @overload
    def __call__(self, /, upd: Updater[T], *, _ui_update: bool = False) -> None:
        ...

    def __call__(
        self, /, upd: T | Updater[T] | Nothing = Nothing.x, *, _ui_update: bool = False
    ) -> T | None:
        assert ctx.in_tx

        if upd is Nothing.x:
            assert ctx.cur_comp is not None

            self._listeners[ctx.cur_comp.id] = ctx.cur_comp
            ctx.cur_comp.signals[self.id] = self

            return self._value

        self._updates.append(upd)
        ctx.updated_signals[self.id] = self
        if not _ui_update:
            ctx.signals_updated_from_code[self.id] = self
        else:
            self._ui_update = True

        self._mark_listeners()

        return None

    def _mark_listeners(self) -> None:
        for x in self._listeners.values():
            x.stale = True
            ctx.stale_nodes[x.id] = x

    def _apply_updates(self) -> None:
        for upd in self._updates:
            if isinstance(upd, Updater):
                self._value = upd.f(self._value)
                return

            self._value = upd

        self._updates = []

    def __repr__(self) -> str:
        return f"{self._name}@{self.id}"


def global_var_signal(key: str) -> Signal[object] | None:
    return _inject.kernel.k_globals.get_signal(key)
