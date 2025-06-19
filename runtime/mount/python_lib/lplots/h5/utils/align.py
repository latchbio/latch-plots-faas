import traceback
from collections.abc import Awaitable, Callable
from contextlib import redirect_stderr, redirect_stdout
from enum import Enum
from io import BytesIO, StringIO
from types import ModuleType

import numpy as np
from lplots.h5.utils import auto_install
from PIL import Image
from scipy.linalg import lstsq

ad = auto_install.ad


class AlignmentMethod(Enum):
    affine = "affine"
    stalign = "stalign"


async def capture_output(blocking_work: Callable[[], None], on_progress:
                         Callable[[object], Awaitable[None]], stage: str,
                         progress_msg_base: dict[str, any]) -> bool:
    error_info: Exception | None = None
    buf_out, buf_err = StringIO(), StringIO()
    try:
        with redirect_stdout(buf_out), redirect_stderr(buf_err):
            blocking_work()
            return False
    except Exception as exc:
        e = exc
        error_info = traceback.format_exc()
        return True
    finally:
        await on_progress(
                {**progress_msg_base,
                 "value": {
                 "stage": stage,
                 "stdout": buf_out.getvalue(),
                 "stderr": buf_err.getvalue(),
                 "error": error_info,
                }})
        if error_info is not None:
            raise e


async def align_image(
        scatter_data_key: str,
        new_scatter_data_key: str,
        points_I: list[list[float]],
        points_J: list[list[float]],
        alignment_method: AlignmentMethod,
        image_bytes: bytes | None,
        adata: ad.AnnData,
        widget_session_key: str,
        on_progress: Callable[[object], Awaitable[None]]
        ) -> None:

    progress_msg_base = {"type": "h5", "op": "align_image", "progress": True,
                         "key": widget_session_key}

    # Each (k, 2)
    # k := # anchor points
    points_I = np.asarray(points_I, dtype=float)
    points_J = np.asarray(points_J, dtype=float)

    k = points_J.shape[0]
    ones = np.ones((k, 1))
    A = np.hstack([points_J, ones])              # (k,3)
    B = points_I                                 # (k,2)

    L: np.ndarray | None = None   # outer scope
    T: np.ndarray | None = None

    # Least-squares:  A · M ≈ B   where M (3, 2)
    def lstsq_work() -> None:
        nonlocal L, T
        M, *_ = lstsq(A, B)
        L = M[:2, :].T      # (2, 2) : linear transformation
        T = M[2, :]         # (2,) : translation

    await capture_output(lstsq_work, on_progress, "lstsq", progress_msg_base)

    X_data = adata.obsm[scatter_data_key]            # (N,2) x,y
    xs, ys = X_data.T

    if alignment_method == AlignmentMethod.affine.value:
        J = np.vstack([xs, ys])             # (2,N)

        I = (L @ J + T[:, None])            # (2,N)

        X_aligned = np.column_stack([I[0], I[1]])  # (N,2)

        adata.obsm[new_scatter_data_key] = X_aligned
        return

    assert image_bytes is not None

    STalign: ModuleType | None = None  # type: ignore
    torch: ModuleType | None = None  # type: ignore

    def install_and_import_work() -> None:
        nonlocal STalign, torch
        STalign, torch = auto_install.install_and_import_stalign_and_torch()

    await capture_output(install_and_import_work, on_progress,
                         "install-and-import", progress_msg_base)

    V = np.array(Image.open(BytesIO(image_bytes)).convert("RGB")) / 255.0
    I_img: np.ndarray | None = None        # (3, H_I, W_I) float32 in [0, 1]

    def normalize_work() -> None:
        nonlocal I_img

        I_img = STalign.normalize(V).transpose(2, 0, 1)   # RGB → CHW

    await capture_output(normalize_work, on_progress, "normalize",
                         progress_msg_base)

    I_y = np.arange(I_img.shape[1])
    I_x = np.arange(I_img.shape[2])

    J_x: np.ndarray | None = None      # (W_J,) grid coordinates
    J_y: np.ndarray | None = None      # (H_J,) grid coordinates
    M:   np.ndarray | None = None      # (1, H_J, W_J) marker image

    def rasterize_work() -> None:
        nonlocal J_x, J_y, M
        J_x, J_y, M, _ = STalign.rasterize(xs, ys, dx=5)

    await capture_output(rasterize_work, on_progress, "rasterize",
                         progress_msg_base)

    # M shape (1, H_J, W_J) but dummy channel
    M2 = M.squeeze(0)

    J_img: np.ndarray | None = None

    def normalize_work() -> None:
        nonlocal J_img
        J_img = STalign.normalize(np.stack([M2] * 3))         # make RGB-like 3×H_J×W_J

    await capture_output(normalize_work, on_progress, "normalize",
                         progress_msg_base)

    if torch.cuda.is_available():
        device = "cuda:0"
    else:
        device = "cpu"

    # https://nbviewer.org/urls/curioanalysisbioinformatics.s3.us-west-1.amazonaws.com/STalign/STalign_Seeker_mouse_embryo_h%26e.ipynb
    params = {"L": L,
              "T": T,
              "niter": 50,
              "pointsI": points_I,
              "pointsJ": points_J,
              "device": device,
              "sigmaM": 0.15,
              "sigmaB": 0.10,
              "sigmaA": 0.11,
              "epV": 10,
              "muB": torch.tensor([0, 0, 0]),
              "muA": torch.tensor([1, 1, 1])
    }

    out: dict[str, any] | None = None

    def lddmm_work() -> None:
        nonlocal out
        out = STalign.LDDMM(
            [I_y, I_x], I_img,
            [J_y, J_x], J_img,
            **params
        )

    await capture_output(lddmm_work, on_progress, "lddmm", progress_msg_base)

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
    adata.obsm[new_scatter_data_key] = tpts.cpu().numpy()
