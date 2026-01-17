# Notebook Cells for UMI Count Filtering, Total cells: 4


## Tab Marker [DEFAULT]
TAB_NAME: Tab 1
TAB_ID: DEFAULT
TYPE: Default Tab Marker
---

## Cell [0] (in Tab 1)
BELONGS_TO_TAB: Tab 1
CELL_ID: cid:0@18199335914055971290:Map
CELL_INDEX: 0
TYPE: code
STATUS: ran
CODE_CELL_ID: 127016
CODE_START
import scanpy as sc
from latch.ldata.path import LPath
from pathlib import Path

# Load the background-removed h5ad file
lp = LPath("latch://38438.account/curio_ovary/v3/ovary_1hr_background_removed.h5ad")
fname = lp.name()
suffix = Path(fname).suffix if fname else ""
local_p = Path(f"{lp.node_id()}{suffix}")
lp.download(local_p, cache=True)

adata = sc.read_h5ad(local_p)
print(f"Loaded data shape: {adata.shape}")
print(f"Beads before filtering: {adata.n_obs}")

# Check UMI distribution
if 'total_counts' not in adata.obs.columns:
    sc.pp.calculate_qc_metrics(adata, inplace=True)

print(f"\nUMI count distribution:")
print(adata.obs['total_counts'].describe())
CODE_END

REACTIVITY:
- Signals defined: adata
- Depends on signals: None
- Depends on cells: None

## Cell [1] (in Tab 1)
BELONGS_TO_TAB: Tab 1
CELL_ID: cid:672@18199335914055971290:Map
CELL_INDEX: 1
TYPE: code
STATUS: ran
CODE_CELL_ID: 127017
CODE_START
from lplots.widgets.text import w_text_output
import json

# Apply minimum UMI count filter
# Standard threshold: minimum 100 UMIs per bead
min_umi_threshold = 100

beads_before = adata.n_obs
adata_filtered = adata[adata.obs['total_counts'] >= min_umi_threshold].copy()
beads_after = adata_filtered.n_obs

# Calculate filtering statistics
beads_removed = beads_before - beads_after
percent_retained = (beads_after / beads_before) * 100

print(f"Minimum UMI threshold: {min_umi_threshold}")
print(f"Beads before filtering: {beads_before}")
print(f"Beads after filtering: {beads_after}")
print(f"Beads removed: {beads_removed}")
print(f"Percent retained: {percent_retained:.1f}%")

# Store result
result = {"beads_after_filtering": beads_after}
print(f"\nResult: {json.dumps(result)}")

w_text_output(
    content=f"✓ Filtered for minimum UMI count ≥ {min_umi_threshold}\n\n**Beads after filtering: {beads_after:,}** ({percent_retained:.1f}% retained)",
    appearance={"message_box": "success"}
)
CODE_END

WIDGETS:
- WIDGET: text_output |  | 127017/0

REACTIVITY:
- Signals defined: adata_filtered, beads_after, beads_before, beads_removed, json, min_umi_threshold, percent_retained, result, w_text_output
- Depends on signals: None
- Depends on cells: None

## Cell [2] (in Tab 1)
BELONGS_TO_TAB: Tab 1
CELL_ID: cid:1671@18199335914055971290:Map
CELL_INDEX: 2
TYPE: code
STATUS: ran
CODE_CELL_ID: 127018
CODE_START
import scanpy as sc
from latch.ldata.path import LPath
from pathlib import Path

# Load the raw Seeker h5ad file
lp = LPath("latch://38438.account/curio_ovary/v3/GSM7689281_adata_ovary_1hr_spatial_raw_counts.h5ad")
fname = lp.name()
suffix = Path(fname).suffix if fname else ""
local_p = Path(f"{lp.node_id()}{suffix}")
lp.download(local_p, cache=True)

adata_raw = sc.read_h5ad(local_p)
print(f"Loaded data shape: {adata_raw.shape}")
print(f"Total beads: {adata_raw.shape[0]}")
print(f"Total genes: {adata_raw.shape[1]}")

# Calculate QC metrics if not present
if 'total_counts' not in adata_raw.obs.columns:
    sc.pp.calculate_qc_metrics(adata_raw, inplace=True)

print(f"\nUMI count distribution:")
print(adata_raw.obs['total_counts'].describe())
CODE_END

REACTIVITY:
- Signals defined: LPath, Path, adata_raw, fname, local_p, lp, sc, suffix
- Depends on signals: None
- Depends on cells: None

## Cell [3] (in Tab 1)
BELONGS_TO_TAB: Tab 1
CELL_ID: cid:3913@18199335914055971290:Map
CELL_INDEX: 3
TYPE: code
STATUS: ran
CODE_CELL_ID: 127019
CODE_START
import sys
sys.path.insert(0, "/opt/latch/plots-faas/runtime/mount/agent_config/context/technology_docs/takara/lib")

from takara import remove_background, KitType
import numpy as np

# First, let's look at the UMI distribution to determine the threshold
import matplotlib.pyplot as plt

log10_umi = np.log10(adata_raw.obs["total_counts"] + 1)
fig, ax = plt.subplots(figsize=(10, 5))
ax.hist(log10_umi, bins=100, edgecolor='black')
ax.set_xlabel("log10(UMI + 1)")
ax.set_ylabel("Count")
ax.set_title("UMI Distribution - Choose threshold at valley between background and tissue")
ax.axvline(x=1.4, color='red', linestyle='--', label='Default threshold (1.4)')
ax.legend()
plt.tight_layout()

from lplots.widgets.plot import w_plot
w_plot(source=fig, label="UMI Distribution for Background Removal")
plt.close()

print("Default min_log10_umi threshold: 1.4")
print("Adjust if needed based on the histogram above")
CODE_END

REACTIVITY:
- Signals defined: sys
- Depends on signals: None
- Depends on cells: None