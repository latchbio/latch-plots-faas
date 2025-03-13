import inspect
import sys
import uuid
from collections.abc import AsyncGenerator, Awaitable, Callable
from contextlib import asynccontextmanager
from dataclasses import dataclass, field
from traceback import print_exc
from typing import Any, Generic, Self, TextIO, TypeAlias, TypeVar, overload

import dill

from . import _inject
from .persistence import SerializedNode, SerializedSignal
from .utils.nothing import Nothing
from .widgets._emit import WidgetState

T = TypeVar("T")
R = TypeVar("R")
Computation: TypeAlias = Callable[..., R]

live_nodes: dict[str, "Node"] = {}
live_node_names: set[str] = set()


def graphviz() -> None:
    from pathlib import Path

    with Path("graph.dot").open("w") as f:
        f.write("digraph G {\n")

        for n in live_nodes.values():
            n.graphviz(f)

        f.write("}\n")


@dataclass(kw_only=True)
class Node:
    f: Computation[Any]
    stale: bool = False
    disposed: bool = False

    parent: Self | None
    children: dict[str, Self] = field(default_factory=dict)

    signals: dict[str, "Signal"] = field(default_factory=dict)

    cell_id: str | None = None
    name: str | None = None
    _id: str

    widget_states: dict[str, WidgetState] = field(default_factory=dict)
    widget_state_idx = 0

    @property
    def id(self) -> str:
        assert self._id is not None
        return self._id

    def __post_init__(self) -> None:
        if self._id is None:
            self._id = str(uuid.uuid4())

        live_nodes[self.id] = self

        if self.parent is not None:
            self.parent.children[self.id] = self
            self.cell_id = self.parent.cell_id
        elif self.name is None:
            self.name = self.cell_id

        if self.name is None:
            self.name = self.f.__name__

        if self.name in live_node_names:
            raise ValueError(f"reactive node name is not unique: {self.name!r}")

        live_node_names.add(self.name)

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
            if self.name is not None:
                live_node_names.remove(self.name)

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
    signals_update_from_code: dict[str, "Signal"] = field(default_factory=dict)
    stale_nodes: dict[str, Node] = field(default_factory=dict)

    in_tx: bool = False

    async def run(
        self, f: Callable[..., Awaitable[R]], *, _cell_id: str | None = None
    ) -> R:
        async with self.transaction:
            self.cur_comp = Node(f=f, parent=self.cur_comp, cell_id=_cell_id, _id=None)

            try:
                if inspect.iscoroutinefunction(f):
                    return await f()
                return f()
            finally:
                self.cur_comp = self.cur_comp.parent

    async def _tick(self) -> None:
        tick_updated_signals = self.signals_update_from_code
        self.signals_update_from_code = {}

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

            async with self.transaction:
                for n, p in to_dispose.values():
                    self.cur_comp = p

                    # we do this here rather than in self.run
                    # because the node is initialized by the cell function
                    # therefor it would not have a cell_id set until before
                    # self.run enters f()
                    #
                    # for this to work with top level nodes, the kernel calls
                    # set_active_cell manually before calling ctx.run()
                    if n.cell_id is not None:
                        await _inject.kernel.set_active_cell(n.cell_id)

                    try:
                        await self.run(n.f, _cell_id=n.cell_id)
                    except Exception:
                        print_exc()
                    finally:
                        self.cur_comp = None
        finally:
            await _inject.kernel.on_tick_finished(tick_updated_signals)

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

    def __init__(
        self, initial: T, *, name: str | None = None, _id: str | None = None
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

    @property
    def id(self) -> str:
        assert self._id is not None
        return self._id

    def sample(self) -> T:
        return self._value

    def serialize(self) -> SerializedSignal:
        error_msg = None
        try:
            val = dill.dumps(self._value)
        except Exception as e:
            val = None
            error_msg = f"Failed to pickle {self._name}: {e}"
        return SerializedSignal(
            value=val,
            name=self._name,
            listeners=self._listeners.keys(),
            error_msg=error_msg,
        )

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
            ctx.signals_update_from_code[self.id] = self

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

    # todo(maximsmol): dispose of signals too to avoid memory leaks

    def __repr__(self) -> str:
        return f"{self._name}@{self.id}"


K = TypeVar("K")
V = TypeVar("V")


class RDict(dict[K, Signal[V]]):
    def __setitem__(self, k: K, v: V) -> None:
        # print(f"SET {k} = {v}")
        if k not in self:
            super().__setitem__(k, Signal(v, name=str(k)))

        super().__getitem__(k)(v)

    def __getitem__(self, k: K) -> V:
        # print(f"GET {k}")
        return super().__getitem__(k)()


def global_var_signal(key: str) -> Signal[object] | None:
    return _inject.kernel.k_globals.get_signal(key)
