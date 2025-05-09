import subprocess  # noqa: S404
import sys
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    import anndata  # type: ignore

if TYPE_CHECKING:
    import STalign  # type: ignore
    import torch


def install_and_import_anndata() -> "anndata":  # type: ignore
    try:
        import anndata as ad  # type: ignore  # noqa: PLC0415
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "anndata"])
        import anndata as ad  # type: ignore  # noqa: PLC0415
    return cast("anndata", ad)


ad = install_and_import_anndata()


def install_and_import_stalign_and_torch() -> tuple["STalign", "torch"]:  # type: ignore
    try:
        import STalign  # type: ignore  # noqa: PLC0415
    except ImportError:

        # note(kenny): author's dumped pip freeze into build requirements in
        # most restrictive way. Only way I was able to install.
        subprocess.check_call([
            sys.executable,
            "-m", "pip", "install", "--no-deps",
            "git+https://github.com/JEFworks-Lab/STalign.git@b2068edc98974efa54537eca194736e177bbe11d"
        ])
        import STalign  # type: ignore  # noqa: PLC0415
        # note(kenny): need strict installation order and dependency resolution
        # might violate that in a single call
        subprocess.check_call([sys.executable, "-m", "pip", "install", "torch===2.7.0"])
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pynrrd==1.1.3"])
        subprocess.check_call([sys.executable, "-m", "pip", "install", "tornado==6.4.2"])
        import torch  # type: ignore  # noqa: PLC0415

    return cast("STalign", STalign), cast("torch", torch)


stalign, torch = install_and_import_stalign_and_torch()
