import asyncio
import contextlib
from collections import OrderedDict, defaultdict
from dataclasses import dataclass

import orjson
from latch_asgi.context.websocket import Context
from latch_asgi.framework.websocket import WebsocketConnectionClosedError


@dataclass(frozen=True, kw_only=True)
class UserProfile:
    auth0_sub: str | None
    connection_idx: int
    name: str
    picture_url: str | None


async def try_send_message(ctx: Context, msg: str) -> None:
    with contextlib.suppress(WebsocketConnectionClosedError):
        await ctx.send_message(msg)


class PlotsContextManager:
    def __init__(self) -> None:
        self.contexts: OrderedDict[str, tuple[Context, UserProfile]] = (
            OrderedDict()
        )  # session_hash -> (context, UserProfile)
        self.session_count_by_user = defaultdict(int)
        self.session_owner: str | int | None = None
        self.unique_users: set[UserProfile] = set()

    async def add_context(
        self,
        sess_hash: str,
        context: Context,
        connection_idx: int,
        auth0_sub: str | None = None,
        picture_url: str | None = None,
        name: str | None = None,
    ) -> None:
        use_auth0 = (
            auth0_sub is not None and picture_url is not None and name is not None
        )
        user_key = auth0_sub if use_auth0 else connection_idx
        assert user_key is not None

        name = name if name is not None else f"Anonymous {connection_idx}"
        user = UserProfile(
            auth0_sub=auth0_sub,
            connection_idx=connection_idx,
            picture_url=picture_url,
            name=name,
        )

        self.contexts[sess_hash] = (context, user)
        self.unique_users.add(user)

        self.session_count_by_user[user_key] += 1
        if self.session_owner is None:
            self.session_owner = user_key

        await self.broadcast_users()

    async def delete_context(self, session_hash: str) -> None:
        if session_hash not in self.contexts:
            return

        _, user = self.contexts[session_hash]
        del self.contexts[session_hash]

        user_key = user.auth0_sub if user.auth0_sub is not None else user.connection_idx
        self.session_count_by_user[user_key] -= 1

        if self.session_count_by_user[user_key] == 0 and self.session_owner == user_key:
            self.unique_users.remove(user)
            self.session_owner = next(iter(self.contexts.values()))[1].auth0_sub

        await self.broadcast_users()

    async def broadcast_message(self, msg: str) -> None:
        tasks = [
            asyncio.create_task(try_send_message(ctx, msg))
            for ctx, _ in self.contexts.values()
        ]
        await asyncio.gather(*tasks)

    async def broadcast_users(self) -> None:
        await self.broadcast_message(
            orjson.dumps({"type": "users", "users": list(self.unique_users)}).decode()
        )
