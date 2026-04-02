from __future__ import annotations

from typing import Any, TypeAlias, TypedDict

import orjson
from latch_data_validation.data_validation import validate
from utils import auth_token_sdk, gql_query, pod_id

Palettes: TypeAlias = dict[str, list[dict[str, Any]]]


class _PlotNotebook(TypedDict):
    id: str
    metadata: str | None


class _PodInfo(TypedDict):
    plotNotebook: _PlotNotebook | None


class _PalettesData(TypedDict):
    podInfo: _PodInfo


class _PalettesResp(TypedDict):
    data: _PalettesData


def _default() -> Palettes:
    return {"categorical": [], "continuous": []}


async def get() -> Palettes:
    if pod_id is None:
        return _default()

    resp = await gql_query(
        query="""
            query GetNotebookPalettes($podId: BigInt!) {
                podInfo(id: $podId) {
                    plotNotebook {
                        id
                        metadata
                    }
                }
            }
        """,
        variables={"podId": pod_id},
        auth=auth_token_sdk,
    )

    data = validate(resp, _PalettesResp)["data"]
    plot_notebook = data["podInfo"]["plotNotebook"]

    if plot_notebook is None:
        return _default()

    metadata_str = plot_notebook["metadata"]
    if metadata_str is None:
        return _default()

    metadata = orjson.loads(metadata_str) if isinstance(metadata_str, str) else metadata_str

    palettes: Palettes = _default()

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
