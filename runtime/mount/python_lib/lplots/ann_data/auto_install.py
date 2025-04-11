import subprocess  # noqa: S404
import sys
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    import anndata  # type: ignore


def install_and_import_anndata() -> "anndata":  # type: ignore
    try:
        import anndata as ad  # type: ignore  # noqa: PLC0415
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "anndata"])
        import anndata as ad  # type: ignore  # noqa: PLC0415
    return cast("anndata", ad)


ad = install_and_import_anndata()
