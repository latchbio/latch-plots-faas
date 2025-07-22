import asyncio
import base64
import json
import sys
import tempfile
import traceback
from collections.abc import Awaitable, Callable
from enum import Enum
from pathlib import Path

import numpy as np

from lplots.h5.utils import auto_install

ad = auto_install.ad


class AlignmentMethod(Enum):
    affine = "affine"
    stalign = "stalign"


async def align_image(
        scatter_data_key: str,
        new_scatter_data_key: str,
        points_i: list[list[float]],
        points_j: list[list[float]],
        alignment_method: str,
        image_bytes: bytes | None,
        adata: ad.AnnData,
        widget_session_key: str,
        send: Callable[[object], Awaitable[None]]
        ) -> None:

    try:
        if scatter_data_key not in adata.obsm:
            await send({
                "type": "h5",
                "op": "align_image",
                "progress": True,
                "key": widget_session_key,
                "value": {
                    "data": {
                        "stage": "validation",
                        "error": f"Scatter data key '{scatter_data_key}' not found in adata.obsm"
                    }
                }
            })
            return

        scatter_data = np.asarray(adata.obsm[scatter_data_key]).tolist()

        subprocess_input = {
            "scatter_data": scatter_data,
            "points_i": points_i,
            "points_j": points_j,
            "alignment_method": alignment_method,
            "widget_session_key": widget_session_key
        }

        if image_bytes is not None:
            subprocess_input["image_bytes_b64"] = base64.b64encode(image_bytes).decode("utf-8")

        subprocess_script = Path(__file__).parent / "align_subprocess.py"

        # temp file to avoid "Argument list too long" error
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=True, encoding="utf-8") as temp_file:
            json.dump(subprocess_input, temp_file)
            temp_file.flush()
            temp_file_path = temp_file.name

            process = await asyncio.create_subprocess_exec(
                sys.executable,
                str(subprocess_script),
                temp_file_path,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )

            aligned_coordinates = None

            async def monitor_stdout() -> None:
                nonlocal aligned_coordinates
                if process.stdout is None:
                    return

                while True:
                    line = await process.stdout.readline()
                    if not line:
                        break

                    try:
                        message = json.loads(line.decode().strip())

                        if message["type"] == "progress":
                            await send({
                                "type": "h5",
                                "op": "align_image",
                                "progress": True,
                                "key": widget_session_key,
                                "value": {
                                    "data": {
                                        "stage": message["stage"]
                                    }
                                }
                            })
                        elif message["type"] == "result":
                            aligned_coordinates = message["aligned_coordinates"]
                        elif message["type"] == "error":
                            await send({
                                "type": "h5",
                                "op": "align_image",
                                "progress": True,
                                "key": widget_session_key,
                                "value": {
                                    "data": {
                                        "stage": "subprocess_error",
                                        "error": message["error"]
                                    }
                                }
                            })

                    except json.JSONDecodeError:
                        continue

            stdout_task = asyncio.create_task(monitor_stdout())

            returncode = await process.wait()
            await stdout_task

            if returncode != 0:
                stderr_output = ""
                if process.stderr:
                    stderr_bytes = await process.stderr.read()
                    stderr_output = stderr_bytes.decode()

                await send({
                    "type": "h5",
                    "op": "align_image",
                    "progress": True,
                    "key": widget_session_key,
                    "value": {
                        "data": {
                            "stage": "subprocess_failed",
                            "error": f"Subprocess failed with return code {returncode}. Stderr: {stderr_output}"
                        }
                    }
                })
                return

            if aligned_coordinates is None:
                await send({
                    "type": "h5",
                    "op": "align_image",
                    "progress": True,
                    "key": widget_session_key,
                    "value": {
                        "data": {
                            "stage": "no_result",
                            "error": "Subprocess completed but no result was received"
                        }
                    }
                })
                return

            adata.obsm[new_scatter_data_key] = np.array(aligned_coordinates, dtype=float)

            await send({
                "type": "h5",
                "op": "align_image",
                "progress": True,
                "key": widget_session_key,
                "value": {
                    "data": {
                        "stage": "completed",
                        "success": True
                    }
                }
            })

    except Exception as e:
        error_msg = f"Alignment error: {e!s}\n{traceback.format_exc()}"
        await send({
            "type": "h5",
            "op": "align_image",
            "progress": True,
            "key": widget_session_key,
            "value": {
                "data": {
                    "stage": "error",
                    "error": error_msg
                }
            }
        })
