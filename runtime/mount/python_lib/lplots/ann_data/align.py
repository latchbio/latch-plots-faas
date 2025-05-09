from enum import Enum
from io import BytesIO

import numpy as np
from lplots.ann_data import auto_install
from PIL import Image
from scipy.linalg import lstsq

ad = auto_install.ad


class AlignmentMethod(Enum):
    affine = "affine"
    stalign = "stalign"


def align_image(
        scatter_data_key: str,
        points_I: list[list[float]],
        points_J: list[list[float]],
        alignment_method: AlignmentMethod,
        image_bytes: bytes | None,
        adata: ad.AnnData,
        ) -> str:

    STalign, torch = auto_install.install_and_import_stalign_and_torch()

    # Each (k, 2)
    # k := # anchor points
    points_I = np.asarray(points_I, dtype=float)
    points_J = np.asarray(points_J, dtype=float)

    ones = np.ones((points_J.shape[0], 1))
    A = np.hstack([points_J, ones])              # (k,3)
    B = points_I                                 # (k,2)

    # Least-squares:  A · M ≈ B   where M (3, 2)
    M, *_ = lstsq(A, B)
    L = M[:2, :].T      # (2, 2) : linear transformation
    T = M[2, :]         # (2,) : translation

    X_data = adata.obsm[scatter_data_key]            # (N,2) x,y
    xs, ys = X_data.T

    if alignment_method == AlignmentMethod.affine:
        J = np.vstack([xs, ys])             # (2,N)

        I = (L @ J + T[:, None])            # (2,N)

        X_aligned = np.column_stack([I[0], I[1]])  # (N,2)

        k = "X_affine_aligned"
        adata.obsm[k] = X_aligned
        return k

    assert image_bytes is not None

    V = np.array(Image.open(BytesIO(image_bytes))) / 255.0
    I_img = STalign.normalize(V).transpose(2, 0, 1)
    I_y = np.arange(I_img.shape[1])
    I_x = np.arange(I_img.shape[2])

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
    k = "X_staligned"
    adata.obsm[k] = tpts.cpu().numpy()
    return k
