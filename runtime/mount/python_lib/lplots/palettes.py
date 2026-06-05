from __future__ import annotations

from collections import defaultdict
from typing import Any, TypeAlias, TypedDict

import orjson
from latch_data_validation.data_validation import validate
from utils import auth_token_sdk, gql_query, pod_id

Palettes: TypeAlias = dict[str, list[dict[str, Any]]]


class PlotNotebook(TypedDict):
    id: str
    palettes: str | None


class PodInfo(TypedDict):
    plotNotebook: PlotNotebook | None


class PalettesData(TypedDict):
    podInfo: PodInfo


class PalettesResp(TypedDict):
    data: PalettesData


default_palette = {"categorical": [], "continuous": []}


async def get() -> Palettes:
    if pod_id is None:
        return default_palette

    resp = await gql_query(
        query="""
            query GetNotebookPalettes($podId: BigInt!) {
                podInfo(id: $podId) {
                    plotNotebook {
                        id
                        palettes
                    }
                }
            }
        """,
        variables={"podId": pod_id},
        auth=auth_token_sdk,
    )

    data = validate(resp, PalettesResp)["data"]
    plot_notebook = data["podInfo"]["plotNotebook"]

    if plot_notebook is None:
        return default_palette

    palettes_str = plot_notebook["palettes"]
    if palettes_str is None:
        return default_palette

    metadata = orjson.loads(palettes_str)

    palettes = defaultdict(list)

    for p in metadata.get("categoricalPalettes") or []:
        palettes["categorical"].append({
            "display_name": p.get("displayName", ""),
            "colors": p.get("colors", []),
        })

    for p in metadata.get("continuousPalettes") or []:
        palettes["continuous"].append({
            "display_name": p.get("displayName", ""),
            "colors": p.get("colors", []),
        })

    return palettes
