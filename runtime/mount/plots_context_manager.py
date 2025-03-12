import asyncio
import contextlib
from collections import defaultdict
from dataclasses import dataclass

import orjson
from latch_asgi.context.websocket import Context
from latch_asgi.framework.websocket import WebsocketConnectionClosedError


@dataclass(frozen=True, kw_only=True)
class UserProfile:
    key: str | int | None
    name: str
    picture_url: str | None


async def try_send_message(ctx: Context, msg: str) -> None:
    with contextlib.suppress(WebsocketConnectionClosedError):
        await ctx.send_message(msg)


class PlotsContextManager:
    def __init__(self) -> None:
        self.contexts: dict[
            str, tuple[Context, UserProfile]
        ] = {}  # session_hash -> (context, UserProfile)
        self.session_count_by_user = defaultdict(lambda: 0)
        self.session_owner: str | int | None = None
        self.unique_users: set[UserProfile] = set()

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
        use_auth0 = (
            auth0_sub is not None and picture_url is not None and name is not None
        )
        user_key = auth0_sub if use_auth0 else f"latch_plots:{connection_idx}"
        assert user_key is not None

        name = name if name is not None else f"Anonymous {connection_idx}"
        user = UserProfile(key=user_key, picture_url=picture_url, name=name)

        self.contexts[sess_hash] = (ctx, user)
        self.unique_users.add(user)

        self.session_count_by_user[user_key] += 1
        if self.session_owner is None:
            self.session_owner = user_key

        await self.broadcast_users()

        return user_key

    async def delete_context(self, session_hash: str) -> None:
        _, user = self.contexts.pop(session_hash, (None, None))
        if user is None:
            return

        self.session_count_by_user[user.key] -= 1

        if self.session_count_by_user[user.key] == 0:
            self.unique_users.remove(user)
            if self.session_owner == user.key:
                self.session_owner = next(iter(self.contexts.values()))[1].key
                del self.session_count_by_user[user.key]

        await self.broadcast_users()

    async def broadcast_message(self, msg: str) -> None:
        async with asyncio.TaskGroup() as tg:
            for ctx, _ in self.contexts.values():
                tg.create_task(try_send_message(ctx, msg))

    async def broadcast_users(self) -> None:
        # todo(rteqs): get rid of extra decode
        await self.broadcast_message(
            orjson.dumps(
                {
                    "type": "users",
                    "users": list(self.unique_users),
                    "session_owner": self.session_owner,
                }
            ).decode()
        )
