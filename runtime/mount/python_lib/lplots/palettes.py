from __future__ import annotations

from typing import Any, TypeAlias

import orjson
from utils import auth_token_sdk, gql_query, pod_id

Palettes: TypeAlias = dict[str, list[dict[str, Any]]]


def _default() -> Palettes:
    return {"categorical": [], "continuous": []}


async def _fetch_notebook_id() -> str | None:
    if pod_id is None:
        return None

    resp = await gql_query(
        query="""
            query GetNotebookId($podId: BigInt!) {
                podInfo(id: $podId) {
                    plotNotebook { id }
                }
            }
        """,
        variables={"podId": pod_id},
        auth=auth_token_sdk,
    )

    return (
        resp.get("data", {})
        .get("podInfo", {})
        .get("plotNotebook", {})
        .get("id")
    )


async def get() -> Palettes:
    notebook_id = await _fetch_notebook_id()
    if notebook_id is None:
        return _default()

    resp = await gql_query(
        query="""
            query GetNotebookPalettes($notebookId: BigInt!) {
                plotNotebookInfo(id: $notebookId) {
                    metadata
                }
            }
        """,
        variables={"notebookId": notebook_id},
        auth=auth_token_sdk,
    )

    metadata_str = (
        resp.get("data", {})
        .get("plotNotebookInfo", {})
        .get("metadata")
    )

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
