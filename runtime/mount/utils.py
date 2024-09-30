import sys
from pathlib import Path
from typing import Any

from aiohttp import ClientSession
from latch_sdk_config.latch import config as latch_config
from yarl import URL

sys.path.append(str(Path(__file__).parent.parent.absolute()))
from config import config

latch_p = Path("/root/.latch")
sdk_token = (latch_p / "token").read_text()
auth_token_sdk = f"Latch-SDK-Token {sdk_token}"

sess: ClientSession | None = None


def get_global_http_sess() -> ClientSession:
    global sess
    if sess is None:
        sess = ClientSession()

    return sess


async def get_presigned_url(path: str) -> str:
    endpoint = latch_config.api.data.get_signed_url

    headers = {"Authorization": auth_token_sdk}
    json_data = {"path": path}

    sess = get_global_http_sess()

    async with sess.post(endpoint, headers=headers, json=json_data) as response:
        res = await response.json()

        if response.status != 200:
            err = res["error"]
            msg = f"failed to fetch presigned url(s) for path {path}"
            if response.status == 400:
                raise ValueError(f"{msg}: download request invalid: {err}")
            if response.status == 401:
                raise RuntimeError(f"authorization token invalid: {err}")
            raise RuntimeError(
                f"{msg} with code {res.status_code}: {res.json()['error']}"
            )

        return res["data"]["url"]


async def gql_query(query: str, variables: dict[str, Any], auth: str) -> Any:
    sess = get_global_http_sess()
    async with sess.post(
        URL(f"https://vacuole.{config.domain}") / "graphql",
        headers={
            "Content-Type": "application/json",
            "User-Agent": "LatchPlotsFaas/0.2.0",
            "Authorization": auth,
        },
        json={"query": query, "variables": variables},
    ) as resp:
        resp.raise_for_status()

        res = await resp.json()
        if "errors" in res:
            raise RuntimeError(f"graphql error: {res}")
        return res
