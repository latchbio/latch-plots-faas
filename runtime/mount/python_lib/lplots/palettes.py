from __future__ import annotations

from collections import defaultdict
from typing import Any, Literal, TypeAlias, TypedDict

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


class Palette(TypedDict):
    colors: list[str]
    type: Literal["categorical", "continuous"]
    display_name: str


latch_palettes: dict[str, Palette] = {
    "basic": {
        "display_name": "Basic",
        "type": "categorical",
        "colors": ["#0634FF", "#FF9301", "#04C700", "#C43BF2", "#FF2700", "#010101"],
    },
    "colorBlindSafe": {
        "display_name": "Color Blind Safe",
        "type": "categorical",
        "colors": ["#000000", "#117F80", "#41007F", "#A966FF", "#FF0066", "#66CBFE"],
    },
    "retroMetro": {
        "display_name": "Retro Metro",
        "type": "categorical",
        "colors": [
            "#ea5545",
            "#f46a9b",
            "#ef9b20",
            "#edbf33",
            "#ede15b",
            "#bdcf32",
            "#87bc45",
            "#27aeef",
            "#b33dc6",
        ],
    },
    "dutchField": {
        "display_name": "Dutch Field",
        "type": "categorical",
        "colors": [
            "#e60049",
            "#0bb4ff",
            "#50e991",
            "#e6d800",
            "#9b19f5",
            "#ffa300",
            "#dc0ab4",
            "#b3d4ff",
            "#00bfa0",
        ],
    },
    "stallion": {
        "display_name": "Stallion",
        "type": "categorical",
        "colors": [
            "#C33530",
            "#282E66",
            "#43884A",
            "#7E2F8A",
            "#E48341",
            "#FAE64D",
            "#8E9ECD",
            "#B570A8",
            "#E0C3DA",
            "#9FD3E2",
            "#96C56C",
            "#E38180",
            "#9584B9",
            "#C25434",
            "#63B9A8",
            "#694D99",
            "#33707A",
            "#731F1C",
            "#D0A970",
            "#3D3D3D",
        ],
    },
    "paired": {
        "display_name": "Paired",
        "type": "categorical",
        "colors": [
            "#a6cee3",
            "#1f78b4",
            "#b2df8a",
            "#33a02c",
            "#fb9a99",
            "#e31a1c",
            "#fdbf6f",
            "#ff7f00",
            "#cab2d6",
            "#6a3d9a",
            "#ffff99",
            "#b15928",
        ],
    },
    "dusky": {
        "display_name": "Dusky Divergent",
        "type": "continuous",
        "colors": [
            "#440154",
            "#482878",
            "#3e4989",
            "#31688e",
            "#26828e",
            "#1f9e89",
            "#35b779",
            "#6ece58",
            "#b5de2b",
            "#fde725",
        ],
    },
    "horizonExtra": {
        "display_name": "Horizon Extra",
        "type": "continuous",
        "colors": [
            "#301437",
            "#572385",
            "#6066b6",
            "#7ca2c2",
            "#c8d0d6",
            "#dacac3",
            "#c6896c",
            "#a84750",
            "#6b1b4d",
            "#2f1436",
        ],
    },
    "sambaNight": {
        "display_name": "Samba Night",
        "type": "continuous",
        "colors": [
            "#000000",
            "#6a009d",
            "#0035dd",
            "#00a4bb",
            "#009b0f",
            "#00e100",
            "#ccf900",
            "#ffb000",
            "#e50000",
            "#cccccc",
        ],
    },
    "solarExtra": {
        "display_name": "Solar Extra",
        "type": "continuous",
        "colors": [
            "#313695",
            "#4a7bb7",
            "#80b7d6",
            "#bde2ee",
            "#eef8df",
            "#feeea5",
            "#fdbf71",
            "#f67b4a",
            "#da372a",
            "#a50026",
        ],
    },
    "whitePurple": {
        "display_name": "White Purple",
        "type": "continuous",
        "colors": [
            "#fcfbfd",
            "#f0eff6",
            "#dfdeed",
            "#c6c7e1",
            "#abaad1",
            "#918dc2",
            "#796eb2",
            "#65489f",
            "#52238d",
            "#3f007d",
        ],
    },
    "whiteBlue": {
        "display_name": "White Blue",
        "type": "continuous",
        "colors": [
            "#f7fbff",
            "#e1edf8",
            "#cbdff1",
            "#abd0e6",
            "#82badb",
            "#59a2cf",
            "#3787c0",
            "#1b6aaf",
            "#084d97",
            "#08306b",
        ],
    },
    "whiteRed": {
        "display_name": "White Red",
        "type": "continuous",
        "colors": [
            "#fff5f0",
            "#fee2d5",
            "#fcc3ac",
            "#fca082",
            "#fb7c5c",
            "#f6553d",
            "#e32f27",
            "#c3161b",
            "#9e0d14",
            "#67000d",
        ],
    },
}

default_categorical_color_palette = latch_palettes["stallion"]["colors"]

default_continuous_color_palette = latch_palettes["solarExtra"]["colors"]
