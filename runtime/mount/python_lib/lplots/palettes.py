from __future__ import annotations

import asyncio
import traceback
from typing import Any

import orjson

from . import _inject

type Palettes = dict[str, list[dict[str, Any]]]


def _default() -> Palettes:
    return {"categorical": [], "continuous": []}


async def _fetch(notebook_id: str | None) -> Palettes:
    if notebook_id is None:
        return _default()

    from utils import auth_token_sdk, gql_query

    try:
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

    except Exception:
        traceback.print_exc()
        return _default()


def get() -> Palettes:
    notebook_id = _inject.kernel.notebook_id
    return asyncio.get_event_loop().run_until_complete(_fetch(notebook_id))
