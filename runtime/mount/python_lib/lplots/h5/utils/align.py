from collections.abc import Awaitable, Callable
from contextlib import redirect_stderr, redirect_stdout
from enum import Enum
from io import BytesIO, StringIO

import numpy as np
from PIL import Image
from scipy.linalg import lstsq

from lplots.h5.utils import auto_install

ad = auto_install.ad


class AlignmentMethod(Enum):
    affine = "affine"
    stalign = "stalign"


progress_msg_base = {"type": "h5", "op": "align_image"}


async def capture_output(blocking_work: Callable[any, None], on_progress:
                         Callable[[object], Awaitable[None]], stage: str) -> None:
    buf_out, buf_err = StringIO(), StringIO()
    try:
        with redirect_stdout(buf_out), redirect_stderr(buf_err):
            blocking_work()
    finally:
        await on_progress(
                {**progress_msg_base,
                 "stage": stage,
                 "stdout": buf_out.getvalue(),
                 "stderr": buf_err.getvalue(),
                })


async def align_image(
        scatter_data_key: str,
        new_scatter_data_key: str,
        points_I: list[list[float]],
        points_J: list[list[float]],
        alignment_method: AlignmentMethod,
        image_bytes: bytes | None,
        adata: ad.AnnData,
        on_progress: Callable[[object], Awaitable[None]]
        ) -> None:

    # Each (k, 2)
    # k := # anchor points
    points_I = np.asarray(points_I, dtype=float)
    points_J = np.asarray(points_J, dtype=float)

    k = points_J.shape[0]
    ones = np.ones((k, 1))
    A = np.hstack([points_J, ones])              # (k,3)
    B = points_I                                 # (k,2)

    # Least-squares:  A · M ≈ B   where M (3, 2)
    def lstsq_work() -> None:
        nonlocal L, T
        M, *_ = lstsq(A, B)
        L = M[:2, :].T      # (2, 2) : linear transformation
        T = M[2, :]         # (2,) : translation

    await capture_output(lstsq_work, on_progress, "lstsq")

    X_data = adata.obsm[scatter_data_key]            # (N,2) x,y
    xs, ys = X_data.T

    if alignment_method == AlignmentMethod.affine.value:
        J = np.vstack([xs, ys])             # (2,N)

        I = (L @ J + T[:, None])            # (2,N)

        X_aligned = np.column_stack([I[0], I[1]])  # (N,2)

        adata.obsm[new_scatter_data_key] = X_aligned
        await on_progress({**progress_msg_base, "stage": "done"})
        return

    assert image_bytes is not None

    await on_progress({**progress_msg_base, "stage": "install-and-import"})

    STalign, torch = auto_install.install_and_import_stalign_and_torch()

    V = np.array(Image.open(BytesIO(image_bytes))) / 255.0
    I_img = STalign.normalize(V).transpose(2, 0, 1)
    I_y = np.arange(I_img.shape[1])
    I_x = np.arange(I_img.shape[2])

    await on_progress({**progress_msg_base, "stage": "rasterize"})

    J_x, J_y, M, _ = STalign.rasterize(xs, ys, dx=5)  # M is (H_J, W_J)

    # M shape (1, H_J, W_J) but dummy channel
    M2 = M.squeeze(0)
    J_img = STalign.normalize(np.stack([M2] * 3))         # make RGB-like 3×H_J×W_J

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

    await on_progress({**progress_msg_base, "stage": "lddm"})

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
    adata.obsm[new_scatter_data_key] = tpts.cpu().numpy()

    await on_progress({**progress_msg_base, "stage": "done"})
