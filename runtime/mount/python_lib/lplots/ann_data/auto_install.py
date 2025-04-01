import subprocess  # noqa: S404
import sys
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    import anndata
    from anndata import AnnData
else:
    anndata = None
    AnnData = None


def install_and_import_anndata() -> "anndata":  # type: ignore  # noqa: PGH003
    try:
        import anndata as ad  # type: ignore  # noqa: PGH003, PLC0415
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "anndata"])
        import anndata as ad  # type: ignore  # noqa: PGH003, PLC0415
    return cast("anndata", ad)


ad = install_and_import_anndata()
