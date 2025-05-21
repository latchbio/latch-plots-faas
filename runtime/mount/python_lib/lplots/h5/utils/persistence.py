import uuid
from pathlib import Path
from typing import TypedDict

from runtime.mount.python_lib.lplots.h5.utils import auto_install

ad = auto_install.ad


class SerializedAnnData(TypedDict):
    key: str
    fname: str
    _is_anndata: bool = True


uns_latch_key_name: str = "_latch_key"


def use_anndata_key(adata: ad.AnnData) -> str:
    key = adata.uns.get(uns_latch_key_name)
    if key is None:
        key = str(uuid.uuid4())
        adata.uns[uns_latch_key_name] = key
    return key


def serialize_anndata(adata: ad.AnnData, snapshot_dir: Path) -> SerializedAnnData:
    key = use_anndata_key(adata)
    fname = f"{key}.h5ad"
    adata_path = snapshot_dir / fname
    if not adata_path.exists():
        adata.write_h5ad(adata_path)
    return SerializedAnnData(key=key, fname=fname)


def load_anndata(s_anndata: SerializedAnnData, snapshot_dir: Path) -> ad.AnnData:
    fname = s_anndata["fname"]
    return ad.read_h5ad(snapshot_dir / fname)
