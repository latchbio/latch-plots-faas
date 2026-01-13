import asyncio
import contextlib
import zlib
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import orjson
from latch_asgi.context.websocket import Context
from latch_asgi.framework.websocket import WebsocketConnectionClosedError

session_count_file = Path("/tmp/user-session-count")  # noqa: S108


@dataclass(frozen=True, kw_only=True)
class UserProfile:
    key: str | None
    name: str
    picture_url: str | None
    is_agent: bool = False


async def try_send_message(ctx: Context, msg: str) -> None:
    with contextlib.suppress(WebsocketConnectionClosedError):
        await ctx.send_message(zlib.compress(msg.encode("utf-8"), level=1))


agent_session_sub = "agent-session"


class PlotsContextManager:
    def __init__(self) -> None:
        self.contexts: dict[
            str, tuple[Context, UserProfile]
        ] = {}  # session_hash -> (context, UserProfile)
        self.session_count_by_user = defaultdict(lambda: 0)
        self.session_owner: str | None = None
        self.unique_users: set[UserProfile] = set()
        self.notebook_id: str | None = None

    def _get_first_non_agent_user_key(self) -> str | None:
        for _, user in self.contexts.values():
            if not user.is_agent:
                return user.key
        return None

    async def add_context(
        self,
        sess_hash: str,
        *,
        ctx: Context,
        connection_idx: int,
        auth0_sub: str | None = None,
        picture_url: str | None = None,
        name: str | None = None,
    ) -> str:
        is_agent = auth0_sub == agent_session_sub

        use_auth0 = (
            auth0_sub is not None and picture_url is not None and name is not None
        )
        user_key = auth0_sub if use_auth0 else f"latch_plots:{connection_idx}"
        assert user_key is not None

        name = name if name is not None else f"Anonymous {connection_idx}"
        user = UserProfile(key=user_key, picture_url=picture_url, name=name, is_agent=is_agent)

        self.contexts[sess_hash] = (ctx, user)
        self.unique_users.add(user)

        self.session_count_by_user[user_key] += 1

        if self.session_owner is None and not is_agent:
            self.session_owner = user_key

        self._write_user_session_count()
        await self.broadcast_users()

        return user_key

    async def delete_context(self, session_hash: str) -> None:
        _, user = self.contexts.pop(session_hash, (None, None))
        if user is None:
            return

        self.session_count_by_user[user.key] -= 1

        if self.session_count_by_user[user.key] == 0:
            self.unique_users.remove(user)
            del self.session_count_by_user[user.key]

            if self.session_owner == user.key:
                if len(self.contexts) == 0:
                    self.session_owner = None
                else:
                    self.session_owner = self._get_first_non_agent_user_key()

        self._write_user_session_count()
        await self.broadcast_users()

    async def broadcast_message(self, msg: str) -> None:
        async with asyncio.TaskGroup() as tg:
            for ctx, _ in self.contexts.values():
                tg.create_task(try_send_message(ctx, msg))

    async def broadcast_users(self) -> None:
        visible_users = [u for u in self.unique_users if not u.is_agent]

        # todo(rteqs): get rid of extra decode
        await self.broadcast_message(
            orjson.dumps({
                "type": "users",
                "users": visible_users,
                "session_owner": self.session_owner,
            }).decode()
        )

    async def override_session_owner(self, user_key: str) -> None:
        self.session_owner = user_key
        await self.broadcast_users()

    def _write_user_session_count(self) -> None:
        user_session_count = sum(
            1 for _, user in self.contexts.values() if not user.is_agent
        )
        session_count_file.parent.mkdir(parents=True, exist_ok=True)
        session_count_file.write_text(str(user_session_count), encoding="utf-8")
