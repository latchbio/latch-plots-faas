# ruff: noqa: N803

import asyncio
import base64
import json
import sys
import tempfile
import traceback
from collections.abc import Awaitable, Callable
from pathlib import Path

import numpy as np

from lplots.h5.utils import auto_install

ad = auto_install.ad


async def align_image(
        scatter_data_key: str,
        new_scatter_data_key: str,
        points_I: list[list[float]],
        points_J: list[list[float]],
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
            "points_I": points_I,
            "points_J": points_J,
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

            process = await asyncio.create_subprocess_exec(
                sys.executable,
                str(subprocess_script.resolve()),
                temp_file.name,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )

            aligned_coordinates = None
            stderr_buffer: str = ""

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

                        msg_type = message.get("type")

                        if msg_type == "progress":
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
                        elif msg_type == "result_file":
                            file_path = message["file_path"]
                            with Path(file_path).open(encoding="utf-8") as fp:  # noqa: ASYNC230
                                aligned_coordinates = json.load(fp)
                            Path(file_path).unlink()
                        elif msg_type == "error":
                            await send({
                                "type": "h5",
                                "op": "align_image",
                                "progress": True,
                                "key": widget_session_key,
                                "value": {
                                    "data": {
                                        "stage": "subprocess",
                                        "error": message["error"]
                                    }
                                }
                            })

                    except json.JSONDecodeError:
                        continue

            async def monitor_stderr() -> None:
                nonlocal stderr_buffer
                if process.stderr is None:
                    return

                while True:
                    line = await process.stderr.readline()
                    if not line:
                        break

                    txt = line.decode(errors="replace")
                    stderr_buffer += txt

                    if len(stderr_buffer) > 5000:
                        stderr_buffer = stderr_buffer[-5000:]

            stdout_task = asyncio.create_task(monitor_stdout())
            stderr_task = asyncio.create_task(monitor_stderr())

            returncode = await process.wait()
            await stdout_task
            await stderr_task

            if returncode != 0:
                await send({
                    "type": "h5",
                    "op": "align_image",
                    "progress": True,
                    "key": widget_session_key,
                    "value": {
                        "data": {
                            "stage": "subprocess",
                            "error": f"Subprocess failed with return code {returncode}. Stderr: {stderr_buffer}"
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
                            "stage": "subprocess",
                            "error": "Subprocess completed but no result was received"
                        }
                    }
                })
                return

            adata.obsm[new_scatter_data_key] = np.array(aligned_coordinates, dtype=float)

            await send({
                "type": "h5",
                "op": "align_image",
                "key": widget_session_key,
                "value": {
                    "data": {
                        "aligned_obsm_key": new_scatter_data_key,
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
