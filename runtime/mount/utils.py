import aiohttp
from latch_cli.utils import get_auth_header
from latch_sdk_config.latch import config


async def get_presigned_url(path: str) -> str:
    endpoint = config.api.data.get_signed_url

    headers = {"Authorization": get_auth_header()}
    json_data = {"path": path}

    async with (
        aiohttp.ClientSession() as session,
        session.post(endpoint, headers=headers, json=json_data) as response,
    ):
        res = await response.json()

        if res.status_code != 200:
            err = res.json()["error"]
            msg = f"failed to fetch presigned url(s) for path {path}"
            if res.status_code == 400:
                raise ValueError(f"{msg}: download request invalid: {err}")
            if res.status_code == 401:
                raise RuntimeError(f"authorization token invalid: {err}")
            raise RuntimeError(
                f"{msg} with code {res.status_code}: {res.json()['error']}"
            )

        data = res.json()

        return data["data"]["url"]
