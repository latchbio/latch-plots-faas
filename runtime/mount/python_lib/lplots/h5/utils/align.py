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

    print(
        f"[align_image] enter: scatter_data_key={scatter_data_key!r} "
        f"new_scatter_data_key={new_scatter_data_key!r} "
        f"alignment_method={alignment_method!r} "
        f"len(points_I)={len(points_I)} len(points_J)={len(points_J)} "
        f"image_bytes_len={len(image_bytes) if image_bytes is not None else None} "
        f"widget_session_key={widget_session_key!r} "
        f"adata.obsm_keys={list(adata.obsm.keys())}",
        file=sys.stdout, flush=True,
    )

    try:
        if scatter_data_key not in adata.obsm:
            print(
                f"[align_image] validation failed: scatter_data_key={scatter_data_key!r} "
                f"not in adata.obsm (keys={list(adata.obsm.keys())})",
                file=sys.stdout, flush=True,
            )
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

        scatter_data_arr = np.asarray(adata.obsm[scatter_data_key])
        scatter_data = scatter_data_arr.tolist()
        print(
            f"[align_image] built scatter_data: shape={scatter_data_arr.shape}",
            file=sys.stdout, flush=True,
        )

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
        print(
            f"[align_image] subprocess_script={subprocess_script} "
            f"exists={subprocess_script.exists()} "
            f"sys.executable={sys.executable}",
            file=sys.stdout, flush=True,
        )

        # temp file to avoid "Argument list too long" error
        # note(aidan): could use socket here. Unclear if it would be better.
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=True, encoding="utf-8") as temp_file:
            json.dump(subprocess_input, temp_file)
            temp_file.flush()
            print(
                f"[align_image] wrote subprocess input to {temp_file.name}",
                file=sys.stdout, flush=True,
            )

            process = await asyncio.create_subprocess_exec(
                sys.executable,
                str(subprocess_script.resolve()),
                temp_file.name,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            print(
                f"[align_image] spawned subprocess: pid={process.pid}",
                file=sys.stdout, flush=True,
            )

            aligned_coordinates = None
            stderr_buffer: str = ""

            async def monitor_stdout() -> None:
                nonlocal aligned_coordinates
                if process.stdout is None:
                    print("[align_image] monitor_stdout: process.stdout is None", file=sys.stdout, flush=True)
                    return

                while True:
                    line = await process.stdout.readline()
                    if not line:
                        print("[align_image] monitor_stdout: stdout EOF", file=sys.stdout, flush=True)
                        break

                    try:
                        message = json.loads(line.decode().strip())

                        msg_type = message.get("type")
                        print(
                            f"[align_image] subprocess stdout msg: type={msg_type!r} "
                            f"keys={list(message.keys())}",
                            file=sys.stdout, flush=True,
                        )

                        if msg_type == "progress":
                            await send({
                                "type": "h5",
                                "op": "align_image",
                                "progress": True,
                                "key": widget_session_key,
                                "value": {
                                    "data": {
                                        "progress_percentage": message["progress_percentage"],
                                        "next_stage": message["next_stage"],
                                        "completed_stage": message.get("completed_stage", None)
                                    }
                                }
                            })
                        elif msg_type == "result_file":
                            file_path = message["file_path"]
                            print(
                                f"[align_image] loading result_file: {file_path}",
                                file=sys.stdout, flush=True,
                            )
                            with Path(file_path).open(encoding="utf-8") as fp:  # noqa: ASYNC230
                                aligned_coordinates = json.load(fp)
                            Path(file_path).unlink()
                            print(
                                f"[align_image] loaded aligned_coordinates: "
                                f"len={len(aligned_coordinates) if aligned_coordinates is not None else None}",
                                file=sys.stdout, flush=True,
                            )
                        elif msg_type == "error":
                            print(
                                f"[align_image] subprocess reported error: {message.get('error')!r}",
                                file=sys.stdout, flush=True,
                            )
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
                        print(
                            f"[align_image] subprocess stdout (non-json): "
                            f"{line.decode(errors='replace').rstrip()!r}",
                            file=sys.stdout, flush=True,
                        )
                        continue

            async def monitor_stderr() -> None:
                nonlocal stderr_buffer
                if process.stderr is None:
                    print("[align_image] monitor_stderr: process.stderr is None", file=sys.stdout, flush=True)
                    return

                while True:
                    line = await process.stderr.readline()
                    if not line:
                        print("[align_image] monitor_stderr: stderr EOF", file=sys.stdout, flush=True)
                        break

                    txt = line.decode(errors="replace")
                    print(
                        f"[align_image] subprocess stderr: {txt.rstrip()!r}",
                        file=sys.stdout, flush=True,
                    )
                    stderr_buffer += txt

                    if len(stderr_buffer) > 5000:
                        stderr_buffer = stderr_buffer[-5000:]

            stdout_task = asyncio.create_task(monitor_stdout())
            stderr_task = asyncio.create_task(monitor_stderr())
            print(
                "[align_image] monitor tasks scheduled; awaiting subprocess exit",
                file=sys.stdout, flush=True,
            )

            returncode = await process.wait()
            print(
                f"[align_image] subprocess exited: returncode={returncode}",
                file=sys.stdout, flush=True,
            )
            await stdout_task
            await stderr_task
            print(
                "[align_image] monitor tasks drained",
                file=sys.stdout, flush=True,
            )

            if returncode != 0:
                print(
                    f"[align_image] non-zero returncode={returncode}; "
                    f"stderr_buffer_tail={stderr_buffer[-1000:]!r}",
                    file=sys.stdout, flush=True,
                )
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
                print(
                    "[align_image] subprocess exited 0 but aligned_coordinates is None",
                    file=sys.stdout, flush=True,
                )
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
            print(
                f"[align_image] wrote adata.obsm[{new_scatter_data_key!r}]; "
                f"shape={adata.obsm[new_scatter_data_key].shape}; sending success",
                file=sys.stdout, flush=True,
            )

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
            print("[align_image] success response sent", file=sys.stdout, flush=True)

    except Exception as e:
        print(
            "[align_image] caught exception:\n" + traceback.format_exc(),
            file=sys.stdout, flush=True,
        )
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
