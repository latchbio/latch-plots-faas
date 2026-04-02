from __future__ import annotations

import os
from pathlib import Path
from typing import Any, TypeAlias

import orjson
from utils import auth_token_sdk, gql_query_sync

Palettes: TypeAlias = dict[str, list[dict[str, Any]]]

_latch_p = Path(os.environ.get("LATCH_SANDBOX_ROOT", "/root/.latch"))
_notebook_id_file = _latch_p / "notebook-id"


def _default() -> Palettes:
    return {"categorical": [], "continuous": []}


def _read_notebook_id() -> str | None:
    try:
        if _notebook_id_file.exists():
            return _notebook_id_file.read_text().strip() or None
    except Exception:
        pass
    return None


def get() -> Palettes:
    notebook_id = _read_notebook_id()
    if notebook_id is None:
        return _default()

    resp = gql_query_sync(
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
