# LPath

**Purpose**

* Use **`LPath`** for all remote `latch://…` paths.
* Use **`pathlib.Path`** **only** for local filesystem.
* Do **not** invent methods beyond those listed here.

**Import**

```python
from latch.ldata.path import LPath
```

## Core Rules

* Remote (`latch://…`) → **LPath** only. Local → **pathlib.Path** only.
* If a widget already returns an **LPath** (e.g., `w_ldata_picker.value`), **use it directly**.
  Don’t wrap `LPath(existing_lpath)`.
* **Always check for `None`** before using any widget value.
* **Always cache file downloads**
* **Never** pass LPath directly to libraries expecting local paths (e.g., `pandas.read_csv`); **download first**.

## Valid Path Forms

* If a path includes a **domain** (e.g., `XXXXX.account`, `XYZ.node`):
  `latch://{domain}/…` (two slashes).
  ❌ `latch:///` with domain (three slashes) is invalid.
* `latch:///…` (three slashes) is only for **relative** (no-domain) paths.

## DO / DON’T

**DO**

```python
# Construct remote
lp = LPath("latch://XXXXX.account/Data/project/file.csv")

# Join remote with "/" operator
child = LPath("latch://XXXXX.account/Data/project") / "file.csv"

# Get remote path string
remote_dir_string = lp.path # returns latch://XXXXX.account/folder 
remote_latch_dir = LatchDir(remote_dir_string)

# Get basename (method, not attribute)
fname = lp.name()

# Download before using local libs
from pathlib import Path
local_p: Path = lp.download()
# e.g., pandas:
import pandas as pd
df = pd.read_csv(local_p)
```

**DON’T**

```python
# ❌ Using Path for remote
from pathlib import Path
bad = Path("latch://…")

# ❌ Re-wrapping an LPath
bad = LPath(lp)
bad = LatchDir(str(lp))
bad = LatchFile(str(lp))

# ❌ Calling str() on an LPath
bad = str(lp)

# ❌ os.path.join or f-strings to build remote paths
bad = os.path.join("latch://…", "file.csv")
bad = f"{lp}/file.csv"

# ❌ Access .name as attribute
bad = lp.name  # use lp.name()
```

## Common Workflow Pattern

```python
# From a widget
from lplots.widgets.ldata import w_ldata_picker
pick = w_ldata_picker(label="Select file/folder")
if pick.value is None:
    exit(0) # or skip safely

lp: LPath = pick.value  # already an LPath

# Inspect and decide how to read
from pathlib import Path

fname = lp.name()
suffix = Path(fname).suffix if fname else ""
local_p = Path(f"{lp.node_id()}{suffix}") 

# Always cache file downloads
lp.download(local_p, cache=True) 

# Example: CSV
import pandas as pd
if fname and fname.endswith(".csv"):
    df = pd.read_csv(local_p)
```

## LPath API (allowed methods & behavior)

* **Metadata / attributes** (lazy; `load_if_missing=True` fetches if needed):

  * `node_id(load_if_missing=True) -> Optional[str]`
  * `name(load_if_missing=True) -> Optional[str]`
  * `type(load_if_missing=True) -> Optional[LDataNodeType]`
  * `size_recursive(load_if_missing=True) -> Optional[int]`
  * `size(load_if_missing=True) -> Optional[int]`
  * `content_type(load_if_missing=True) -> Optional[str]`
  * `is_dir(load_if_missing=True) -> bool`
  * `fetch_metadata() -> None` (forces network fetch & caches)

* **Directory ops**

  * `iterdir() -> Iterator[LPath]` *(network; non-recursive; dir only)*
  * `mkdirp() -> None`
  * `rmr() -> None` *(recursive delete; raises if path missing; network)*

* **Data movement**

  * `copy_to(dst: LPath) -> None` *(remote → remote)*
  * `upload_from(src: pathlib.Path, show_progress_bar: bool=False) -> None`
  * `download(dst: Optional[pathlib.Path]=None, show_progress_bar: bool=False) -> pathlib.Path`

* **Caching downloads (strongly recommended)**

```python
from latch.ldata.path import LPath
from pathlib import Path

lp = LPath("latch:///welcome/deseq2/design.csv")
dst = Path(f"{lp.node_id()}.csv")  # Local cache path
lp.download(dst, cache=True)       # Cached download
```

* **Joining**

  * `LPath / "child.ext"` to create a child remote path.

* **Types**

  * `LDataNodeType`: `account_root | dir | obj | mount | link | mount_gcp | mount_azure`

## Quick Diagnostics

* Getting errors with pandas / seaborn / etc. on `latch://…`?
  → You forgot to **download** to a local `Path` first.
* Seeing `ValueError` from `LPath(...)`?
  → You likely passed an **LPath** into `LPath(...)`. Use it directly.
* Path fails with domain and `latch:///…`?
  → Use **two** slashes with domain: `latch://{domain}/…`.
