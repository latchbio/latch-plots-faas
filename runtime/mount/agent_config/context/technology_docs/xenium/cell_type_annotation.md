# Cell Type Annotation (CellGuide-Based)

This step assigns cell types to clusters by:
1. Using precomputed **cluster-level marker genes** (from DGE).
2. Matching markers to **CellGuide**.
3. Writing the final labels back into `adata.obs`.

---

## 1. High-Level Flow

1. Pick the **cluster column** in `adata.obs` and the matching **DGE key** in `adata.uns`.
2. Convert the DGE result to a tidy DataFrame of markers per cluster.
3. For each cluster, take the top markers and query CellGuide.
4. Aggregate evidence across markers and pick the most supported cell type(s).
5. Store cell-type labels in a new column in `adata.obs`.

---

## 2. CellGuide Databases

CellGuide provides curated mappings between genes, cell types, tissues, and organisms. We use:

- A **per-gene** database to ask “which cell types does this gene mark?”
- A **per-cell-type** database to ask “which genes mark this cell type?”
- A **tissue metadata** file to constrain queries to the correct organism/tissue.

### 2.1 Per-Gene Database

**File:** `cellguide_marker_gene_database_per_gene.json`
Answers: *“Which cell types is this gene a marker for?”*

```python
marker_db_per_gene["Gad1"]["Mus musculus"]["brain"]
```

- Used for automated annotation (main engine of this step).
- Cell types are ordered by marker strength per gene.

### 2.2 Per-Cell-Type Database

**File:** `cellguide_marker_gene_database_per_celltype.json`
Answers: *“What are the marker genes for this cell type?”*

```python
marker_db_per_celltype["Mus musculus"]["brain"]["inhibitory neuron"]
```

- Used for validation, manual review, and marker rationales.

### 2.3 Tissue Options

**File:** `tissues_per_organism.json`
Answers: *“What tissues are available for this organism?”*

```python
tissues_per_organism["Mus musculus"]
```

- Used for validating organism/tissue choices and populating UI dropdowns.

### 2.4 Loading from Latch

```python
import json
from latch.ldata.path import LPath
from pathlib import Path

def load_json_lpath(uri: str):
    lp = LPath(uri)
    local = Path(f"{lp.node_id()}.json")
    lp.download(local, cache=True)
    with local.open() as f:
        return json.load(f)

marker_db_per_gene = load_json_lpath(
    "latch:///cellguide_marker_gene_database_per_gene.json"
)
marker_db_per_celltype = load_json_lpath(
    "latch:///cellguide_marker_gene_database_per_celltype.json"
)
tissues_per_organism = load_json_lpath(
    "latch:///tissues_per_organism.json"
)
```

---

## 2.5 Central Cell-Type Vocabulary Configs

We keep **one vocab config per organism+tissue(+panel)** and a **central index** that points to them.

- The **central index** is what the agent / notebook loads first.
- The **per-tissue vocab configs** (for example brain and kidney) live as separate JSONs in the same docs directory.

A typical layout:

```text
/opt/latch/plots-faas/runtime/mount/agent_config/context/technology_docs/
  xenium/
    cell_type_annotation.md
    cell_type_vocab_index.json
    cell_type_vocab_mus_musculus_brain.json
    cell_type_vocab_mus_musculus_kidney.json
```

### 2.5.1 Central index structure

The central index file `cell_type_vocab_index.json` lists all available vocab configs:

```json
{
  "configs": [
    {
      "organism": "Mus musculus",
      "tissue": "brain",
      "panel_name": "xenium_full_brain",
      "path": "cell_type_vocab_mus_musculus_brain.json"
    }
  ]
}
```

Each entry:

- `organism`: e.g. `"Mus musculus"`
- `tissue`: e.g. `"brain"`, `"kidney"`
- `panel_name`: panel or dataset name (optional)
- `path`: relative filename of the per-tissue vocab JSON in the docs directory

> `panel_name` is optional and primarily used by the agent to disambiguate between multiple vocab configs for the same organism and tissue. Users are not expected to know or set this.

### 2.5.2 Loading the central index and vocab configs (local files)

```python
import json
from pathlib import Path

# Relative path from context directory
DOCS_DIR = Path(__file__).resolve().parent
VOCAB_INDEX_PATH = DOCS_DIR / "cell_type_vocab_index.json"


def load_vocab_index() -> dict:
    with VOCAB_INDEX_PATH.open() as f:
        return json.load(f)


def load_cell_type_vocab_config(
    organism: str,
    tissue: str,
    panel_name: str | None = None,
) -> dict:
    """
    Look up the vocab config for a given (organism, tissue[, panel_name])
    using the central index, then load and return that JSON.

    If panel_name is None, use the first config that matches (organism, tissue).
    """
    index = load_vocab_index()

    # First try exact triple match, if panel_name is provided.
    if panel_name is not None:
        for cfg in index.get("configs", []):
            if (
                cfg.get("organism") == organism
                and cfg.get("tissue") == tissue
                and cfg.get("panel_name") == panel_name
            ):
                vocab_path = DOCS_DIR / cfg["path"]
                with vocab_path.open() as f:
                    return json.load(f)

    # Fallback: any config matching organism + tissue.
    for cfg in index.get("configs", []):
        if (
            cfg.get("organism") == organism
            and cfg.get("tissue") == tissue
        ):
            vocab_path = DOCS_DIR / cfg["path"]
            with vocab_path.open() as f:
                return json.load(f)

    raise ValueError(
        f"No vocab config found for organism={organism}, tissue={tissue}, "
        f"panel_name={panel_name!r}"
    )

```

### 2.5.3 Usage in the annotation flow

After summarizing clusters with CellGuide, load the appropriate vocab config:

```python
vocab_config = load_cell_type_vocab_config(
    organism="Mus musculus",
    tissue="brain",
)
```

The returned `vocab_config` contains:

- `allowed_vocab`: controlled labels for that panel
- `mapping_rules`: regex rules to map raw CellGuide names into the controlled vocabulary

These fields are then used downstream to convert raw CellGuide labels into the final `adata.obs["cell_type"]` labels.

---

## 3. Prerequisites

Your `AnnData` must include:

- **Cluster assignments (required)**
  A column in `adata.obs` whose name follows
  `clustering_{clustering_algorithm}_{clustering_resolution}`
  Examples:
  - `clustering_leiden_0.4`
  - `clustering_louvain_1.0`

- **Expression matrix**
  Expression values available in `adata.X` or an appropriate `adata.layers[...]`
  (already used for DGE upstream).

Additional inputs:

- **CellGuide per-gene database**
  `db[gene][organism][tissue] -> list of candidate cell-type names`
- **Organism and tissue**, for example:
  - `organism = "Mus musculus"`
  - `tissue   = "brain"`

> Tip: Always inspect `adata.obs.columns` and `adata.uns.keys()` before assuming cluster or DGE key names.

---

## 4. Select Clustering and Retrieve Marker Genes

1. Choose a clustering column in `adata.obs` whose name starts with `clustering_`.
2. Choose the matching DGE key in `adata.uns` whose name starts with `rank_genes_groups_clustering_`.
3. Extract marker genes per cluster into a DataFrame.

```python
import scanpy as sc

print("obs columns:", adata.obs.columns.tolist())
print("uns keys:", list(adata.uns.keys()))

cluster_col = "clustering_leiden_0.4"                  # example
dge_key = "rank_genes_groups_clustering_leiden_0.4"    # example

ranked_df = sc.get.rank_genes_groups_df(
    adata,
    group=None,   # all clusters
    key=dge_key,
)
```

---

## 5. Lookup Cell Types in CellGuide

For a list of genes, return CellGuide cell types for the given organism/tissue.

This step operates in **raw CellGuide label space** only. The final labels may
later be normalized into a controlled vocabulary using the vocab configs
described in Section 2.5 (if available).

```python
import json
from pathlib import Path

def lookup_cellguide_celltypes(genes, organism, tissue, db_path):
    db_path = Path(db_path)
    with db_path.open() as f:
        db = json.load(f)

    results = {}
    for g in genes:
        org_entry = db.get(g, {}).get(organism, {})
        cts = org_entry.get(tissue) or org_entry.get("All Tissues")
        if cts:
            results[g] = cts
    return results
```

## 6. Summarize Cluster Annotations

Use the top markers per cluster to vote on the most likely cell type(s).

```python
from collections import Counter
import json
from pathlib import Path

def summarize_clusters(
    ranked_df,
    db_path,
    organism="Mus musculus",
    tissue="brain",
    min_markers=3,
    n_core=10,
):
    db_path = Path(db_path)
    with db_path.open() as f:
        db = json.load(f)

    def _lookup_cellguide_celltypes(genes):
        results = {}
        for g in genes:
            org_entry = db.get(g, {}).get(organism, {})
            cts = org_entry.get(tissue) or org_entry.get("All Tissues")
            if cts:
                results[g] = cts
        return results

    summary = {}

    # Iterate over clusters that actually have DGE results
    for c in sorted(ranked_df["group"].unique()):
        g = ranked_df[ranked_df["group"] == c]
        if g.empty:
            summary[c] = {
                "core_markers": [],
                "most_common_cell_type": [],
                "cell_type_counts": {},
                "markers_for_most_common_cell_type": [],
            }
            continue

        # Use top n_core markers per cluster
        g = g.sort_values("rank")
        core_markers = g["names"].head(n_core).tolist()

        lookup = _lookup_cellguide_celltypes(core_markers)

        # If none of the core markers exist in CellGuide
        if not lookup:
            summary[c] = {
                "core_markers": core_markers,
                "most_common_cell_type": [],
                "cell_type_counts": {},
                "markers_for_most_common_cell_type": [],
            }
            continue

        # Count how many core markers support each cell type
        counter = Counter(ct for v in lookup.values() for ct in v)

        # Only consider cell types supported by at least min_markers core markers
        eligible = {ct: n for ct, n in counter.items() if n >= min_markers}
        if not eligible:
            summary[c] = {
                "core_markers": core_markers,
                "most_common_cell_type": [],
                "cell_type_counts": dict(counter),  # keep raw counts for debugging
                "markers_for_most_common_cell_type": [],
            }
            continue

        max_support = max(eligible.values())
        top_cts = [ct for ct, n in eligible.items() if n == max_support]

        # Markers that support at least one of the top cell types
        top_markers = [
            gene
            for gene, cts in lookup.items()
            if any(ct in top_cts for ct in cts)
        ]

        summary[c] = {
            "core_markers": core_markers,
            "most_common_cell_type": top_cts,
            "cell_type_counts": dict(counter),
            "markers_for_most_common_cell_type": top_markers,
        }

    return summary
```

**Per-cluster outputs:**

- `core_markers`: Top marker genes used for CellGuide matching.
- `most_common_cell_type`: Predicted identity as a list to allow ties or mixed signatures.
- `cell_type_counts`: Number of core markers that support each candidate cell type.
- `markers_for_most_common_cell_type`: Subset of core markers that support at least one top cell type.

These raw CellGuide labels are the **input** to the final labeling step described in cell_type_annotation.md (mapping to `allowed_vocab` when a vocab config is available, or using cleaned CellGuide names otherwise).

---

## Best Practices and Troubleshooting

| Issue | Likely Cause | Recommendation |
|------|--------------|----------------|
| No CellGuide matches | Tissue or organism mismatch | Verify correct organism/tissue keys |
| Missing cluster column | Dataset not clustered | Run UMAP + Leiden/Louvain before annotation |
