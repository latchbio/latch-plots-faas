import functools
import os
import sys
import time
from collections.abc import Awaitable, Callable
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from typing import Any, TypeVar

# Enable by setting LPLOTS_H5_PROFILE=1 in the environment on the dev box.
# Output lines are prefixed with "[h5prof]" so they are easy to grep in logs.
PROFILE_ENABLED = os.environ.get("LPLOTS_H5_PROFILE", "") not in ("", "0", "false")


def _emit(line: str) -> None:
    """Write a profiling line to the process's real stderr (→ systemd journal).

    The kernel reassigns ``sys.stdout``/``sys.stderr`` to a SocketWriter that
    pipes output to the frontend notebook, so a plain ``print`` would NOT reach
    journald. ``sys.__stderr__`` still references the original fd, which systemd
    captures into the journal. Fall back to a raw write on fd 2.
    """
    payload = line if line.endswith("\n") else line + "\n"
    try:
        real_err = sys.__stderr__
        if real_err is not None:
            real_err.write(payload)
            real_err.flush()
            return
    except Exception:
        pass

    try:
        os.write(2, payload.encode("utf-8", errors="replace"))
    except OSError:
        pass


# Unconditional one-time marker at import so we can confirm the new code is
# actually deployed and see whether the flag is on. Remove once debugging done.
_emit(
    f"[h5prof] profiling module loaded "
    f"enabled={PROFILE_ENABLED} "
    f"LPLOTS_H5_PROFILE={os.environ.get('LPLOTS_H5_PROFILE')!r} "
    f"pid={os.getpid()}"
)


@dataclass
class ProfileRun:
    label: str
    started_at: float = field(default_factory=time.perf_counter)
    # ordered: name -> (total_seconds, call_count)
    sections: "OrderedTimings" = field(default_factory=lambda: OrderedTimings())


class OrderedTimings(dict):
    def add(self, name: str, dt: float) -> None:
        total, count = self.get(name, (0.0, 0))
        self[name] = (total + dt, count + 1)


_current_run: ContextVar[ProfileRun | None] = ContextVar(
    "h5_profile_run", default=None
)


@contextmanager
def profile_request(label: str):
    """Top-level span for a single h5 request. Prints a summary on exit."""
    if not PROFILE_ENABLED:
        yield
        return

    run = ProfileRun(label=label)
    token = _current_run.set(run)
    try:
        yield
    finally:
        _current_run.reset(token)
        total_ms = (time.perf_counter() - run.started_at) * 1000
        parts = []
        for name, (secs, count) in run.sections.items():
            ms = secs * 1000
            if count > 1:
                parts.append(f"{name}={ms:.1f}ms(x{count})")
            else:
                parts.append(f"{name}={ms:.1f}ms")
        detail = " ".join(parts)
        _emit(f"[h5prof] {label} total={total_ms:.1f}ms {detail}")


F = TypeVar("F", bound=Callable[..., Awaitable[Any]])


def profile_request_by_op(fn: F) -> F:
    """Decorator for an async handler: opens a profile run labeled `op=<msg op>`.

    Assumes the first positional arg is the message dict (or has `.get("op")`).
    """
    if not PROFILE_ENABLED:
        return fn

    @functools.wraps(fn)
    async def wrapper(*args: Any, **kwargs: Any) -> Any:
        msg = args[0] if args else kwargs.get("msg")
        op = "unknown"
        if isinstance(msg, dict):
            op = str(msg.get("op", "unknown"))

        run = ProfileRun(label=f"op={op}")
        token = _current_run.set(run)
        try:
            return await fn(*args, **kwargs)
        finally:
            _current_run.reset(token)
            total_ms = (time.perf_counter() - run.started_at) * 1000
            parts = []
            for name, (secs, count) in run.sections.items():
                ms = secs * 1000
                if count > 1:
                    parts.append(f"{name}={ms:.1f}ms(x{count})")
                else:
                    parts.append(f"{name}={ms:.1f}ms")
            detail = " ".join(parts)
            _emit(f"[h5prof] {run.label} total={total_ms:.1f}ms {detail}")

    return wrapper  # type: ignore[return-value]


@contextmanager
def measure(label: str):
    """Standalone span that emits its own line immediately on exit.

    Unlike `profile`, this does not require an active `profile_request`, so it
    can be used to time code that runs *after* the request handler returns
    (e.g. response serialization + socket send).
    """
    if not PROFILE_ENABLED:
        yield
        return

    start = time.perf_counter()
    try:
        yield
    finally:
        _emit(f"[h5prof] {label}={(time.perf_counter() - start) * 1000:.1f}ms")


@contextmanager
def profile(name: str):
    """Named sub-span. No-op unless inside an enabled profile_request."""
    run = _current_run.get()
    if run is None:
        yield
        return

    start = time.perf_counter()
    try:
        yield
    finally:
        run.sections.add(name, time.perf_counter() - start)
