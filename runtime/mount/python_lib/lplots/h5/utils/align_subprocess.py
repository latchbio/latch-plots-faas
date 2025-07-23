# ruff: noqa: N803, N806

import base64
import json
import sys
import tempfile
import traceback
from enum import Enum
from io import BytesIO
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.linalg import lstsq

# add mount directory to Python path so we can import lplots
mount_dir = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(mount_dir))

from lplots.h5.utils import auto_install


class AlignmentMethod(Enum):
    affine = "affine"
    stalign = "stalign"


def send_progress_update(completed_stage: str, next_stage: str, progress: float, widget_session_key: str = "") -> None:
    message = {
        "type": "progress",
        "completed_stage": completed_stage,
        "next_stage": next_stage,
        "progress": progress,
        "widget_session_key": widget_session_key
    }
    print(json.dumps(message), file=sys.stdout, flush=True)


def send_result_file(file_path: str, widget_session_key: str = "") -> None:
    message = {
        "type": "result_file",
        "file_path": file_path,
        "widget_session_key": widget_session_key,
    }
    print(json.dumps(message), file=sys.stdout, flush=True)


def send_error(error_message: str, widget_session_key: str = "") -> None:
    message = {
        "type": "error",
        "error": error_message,
        "widget_session_key": widget_session_key
    }
    print(json.dumps(message), file=sys.stdout, flush=True)


def align_image_subprocess(
    scatter_data: list[list[float]],
    points_I: list[list[float]],
    points_J: list[list[float]],
    alignment_method: str,
    image_bytes: bytes | None,
    widget_session_key: str
) -> list[list[float]]:
    points_I_arr = np.asarray(points_I, dtype=float)
    points_J_arr = np.asarray(points_J, dtype=float)
    X_data = np.asarray(scatter_data, dtype=float)  # (N,2) x,y

    send_progress_update("lstsq", "install-and-import", 0.1, widget_session_key)

    k = points_J_arr.shape[0]
    ones = np.ones((k, 1))
    A = np.hstack([points_J_arr, ones])              # (k,3)
    B = points_I_arr                                 # (k,2)

    # Least-squares:  A · M ≈ B   where M (3, 2)
    M, *_ = lstsq(A, B)
    L = M[:2, :].T      # (2, 2) : linear transformation
    T = M[2, :]         # (2,) : translation

    xs, ys = X_data.T

    if alignment_method == AlignmentMethod.affine.value:
        J = np.vstack([xs, ys])             # (2,N)
        I = (L @ J + T[:, None])            # (2,N)
        X_aligned = np.column_stack([I[0], I[1]])  # (N,2)
        return X_aligned.tolist()  # type: ignore[return-value]

    # STalign method
    if image_bytes is None:
        raise ValueError("Image bytes required for STalign method")

    send_progress_update("install-and-import", "normalize_image", 0.2, widget_session_key)
    STalign, torch = auto_install.install_and_import_stalign_and_torch()

    send_progress_update("normalize_image", "rasterize", 0.3, widget_session_key)
    V = np.array(Image.open(BytesIO(image_bytes)).convert("RGB")) / 255.0
    I_img = STalign.normalize(V).transpose(2, 0, 1)   # RGB → CHW

    I_y = np.arange(I_img.shape[1])
    I_x = np.arange(I_img.shape[2])

    send_progress_update("rasterize", "normalize_points", 0.4, widget_session_key)
    J_x, J_y, M, _ = STalign.rasterize(xs, ys, dx=5)

    # M shape (1, H_J, W_J) but dummy channel
    M2 = M.squeeze(0)

    send_progress_update("normalize_points", "lddmm", 0.5, widget_session_key)
    J_img = STalign.normalize(np.stack([M2] * 3))         # make RGB-like 3xH_JxW_J

    if torch.cuda.is_available():
        device = "cuda:0"
    else:
        device = "cpu"

    # https://nbviewer.org/urls/curioanalysisbioinformatics.s3.us-west-1.amazonaws.com/STalign/STalign_Seeker_mouse_embryo_h%26e.ipynb
    params = {"L": L,
              "T": T,
              "niter": 50,
              "pointsI": points_I_arr,
              "pointsJ": points_J_arr,
              "device": device,
              "sigmaM": 0.15,
              "sigmaB": 0.10,
              "sigmaA": 0.11,
              "epV": 10,
              "muB": torch.tensor([0, 0, 0]),
              "muA": torch.tensor([1, 1, 1])
    }

    send_progress_update("lddmm", "transform_points", 0.6, widget_session_key)
    out = STalign.LDDMM(
        [I_y, I_x], I_img,
        [J_y, J_x], J_img,
        **params
    )

    pts = np.stack([xs, ys], axis=1)  # (N,2) in J-coords
    dtype = out["A"].dtype
    pts_torch = torch.tensor(
        pts,
        device=device,
        dtype=dtype
    )

    tpts = STalign.transform_points_source_to_target(
        out["xv"], out["v"], out["A"], pts_torch
    )

    return tpts.cpu().numpy().tolist()


def main() -> None:
    if len(sys.argv) != 2:
        send_error("Usage: align_subprocess.py <json_input_file>")
        sys.exit(1)

    try:
        input_file_path = sys.argv[1]
        with Path(input_file_path).open(encoding="utf-8") as f:
            input_data = json.load(f)

        scatter_data = input_data["scatter_data"]
        points_I = input_data["points_I"]
        points_J = input_data["points_J"]
        alignment_method = input_data["alignment_method"]
        widget_session_key = input_data["widget_session_key"]

        image_bytes = None
        if "image_bytes_b64" in input_data:
            image_bytes = base64.b64decode(input_data["image_bytes_b64"])

        aligned_coordinates = align_image_subprocess(
            scatter_data=scatter_data,
            points_I=points_I,
            points_J=points_J,
            alignment_method=alignment_method,
            image_bytes=image_bytes,
            widget_session_key=widget_session_key,
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False, prefix="aligned_coords_", encoding="utf-8") as tmp:
            json.dump(aligned_coordinates, tmp)
            result_path = tmp.name

        send_result_file(result_path, widget_session_key)

    except Exception as e:
        error_msg = f"Subprocess error: {e!s}\n{traceback.format_exc()}"
        send_error(error_msg, input_data.get("widget_session_key", "") if "input_data" in locals() else "")
        sys.exit(1)


if __name__ == "__main__":
    main()
