import base64
import io
import os
import sys
from enum import Enum
from pathlib import Path
from typing import Any, Literal, NotRequired, TypedDict

from aiohttp import ClientSession
from latch.types.directory import LatchDir
from latch.types.file import LatchFile
from lplots.utils.nothing import _Nothing
from matplotlib.figure import Figure, SubFigure
from typing_extensions import _AnnotatedAlias
from yarl import URL

# todo(rteqs): get rid of this
sys.path.append(str(Path(__file__).parent.parent.absolute()))
from config import config

latch_p = Path(os.environ.get("LATCH_SANDBOX_ROOT", "/root/.latch"))
sdk_token_path = latch_p / "token"
if sdk_token_path.exists():
    sdk_token = sdk_token_path.read_text()
else:
    sdk_token = ""
auth_token_sdk = f"Latch-SDK-Token {sdk_token}"


nucleus_url_path = latch_p / "nucleus-url"
if nucleus_url_path.exists():
    nucleus_url = nucleus_url_path.read_text()
else:
    nucleus_url = "https://nucleus.latch.bio"

sess: ClientSession | None = None

KernelSnapshotStatus = Literal["done", "progress", "start", "error"]


class Trace(TypedDict):
    type: str
    x: str
    y: str
    color_by: NotRequired[str | None]
    error_bar: NotRequired[str | None]
    marker_size: NotRequired[str | None]


class PlotConfig(TypedDict):
    traces: list[Trace]
    custom_data: NotRequired[list[str]]

    facet: NotRequired[str | None]

    xrange: NotRequired[tuple[float | int, float | int]]
    yrange: NotRequired[tuple[float | int, float | int]]

    height_px: NotRequired[float]
    width_px: NotRequired[float]


def get_global_http_sess() -> ClientSession:
    global sess
    if sess is None:
        sess = ClientSession()

    return sess


async def get_presigned_url(path: str) -> str:
    endpoint = "ldata/get-signed-url"

    headers = {"Authorization": auth_token_sdk}
    json_data = {"path": path}

    sess = get_global_http_sess()

    async with sess.post(
        URL(nucleus_url) / endpoint, headers=headers, json=json_data
    ) as response:
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


def plot_to_webp_string(
    fig: Figure | SubFigure, *, quality: int = 80, lossless: bool = False
) -> str:
    buf = io.BytesIO()
    fig.canvas.print_figure(
        buf,
        format="webp",
        dpi=300,
        bbox_inches="tight",
        pil_kwargs={"quality": quality, "lossless": lossless},
    )
    fig_b64 = base64.b64encode(buf.getvalue()).decode("utf-8")
    return f"data:image/webp;base64,{fig_b64}"


def orjson_encoder(obj: object) -> str | None:
    if isinstance(obj, Figure):
        return plot_to_webp_string(obj)

    if hasattr(obj, "figure") and isinstance(obj.figure, Figure):
        fig = obj.figure
        if fig is None:
            return None

        return plot_to_webp_string(fig)

    if isinstance(obj, (LatchFile, LatchDir)):
        return obj.remote_path

    if isinstance(obj, Enum):
        return obj.name

    # todo(rteqs): figure out what's going in orjson interanally when serializing pandas dataframes with Nothign.x
    if isinstance(obj, _Nothing | type | _AnnotatedAlias):
        return str(obj)

    raise TypeError(f"Type {type(obj)} not serializable")
