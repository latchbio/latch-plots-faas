# Marker-based Cell Type Annotation

## Overview
This guide provides a step-by-step workflow for annotating cell types in spatial biology datasets using **marker gene analysis** and the **CellGuide** database.

Two complementary strategies are supported depending on the dataset quality and clustering resolution:

- **Approach 1 — Gene set scoring**  
  Assigns cell types based on expression of predefined marker sets.  
  _Does not require clusters; useful for quick validation or noisy data._

- **Approach 2 — Cluster-based annotation**  
  Identifies differentially expressed genes per cluster and matches them to known cell types.  
  _Requires clear, biologically meaningful clusters._

---

## Critical Decision Point (Agent MUST Ask Before Continuing)

Before any analysis, the agent **MUST ask the user two questions:**

1. “Are you confident in your current clustering quality?”  
2. “Do you prefer a fast, exploratory annotation (gene set scoring) or a slower, more detailed annotation (differential expression + marker matching)?  
   – Fast approach: compute gene-set scores using user-provided markers or CellGuide markers.  
   – Slow approach: run 1-vs-all Wilcoxon differential expression, rank markers, and match them to reference databases.”

Use the user’s answers to select the approach:

- **IF** the user is **not confident** in clustering quality **AND/OR** wants **fast, exploratory annotation**  
  → **THEN run Approach 1 (Gene set scoring)** → Follow documentation in `technology_docs/cell_type_annotation/gene_set_scoring.md`
  - Fast (seconds to minutes), works even with weak clusters, provides coarse preliminary labels. 

- **IF** the user is **confident** in clustering quality **AND** is willing to wait **15–30+ minutes**  
  → **THEN run Approach 2 (Cluster-based annotation)** → Follow documentation in `technology_docs/cell_type_annotation/cluster_based.md`
  - Slow, requires strong clusters. 

The agent MUST NOT proceed until both questions are answered and an approach is explicitly chosen.

---

## Input Data

In AtlasXomics’ processed datasets, the file ending with `_sm_he.h5ad` typically contains gene-activity scores and is the recommended input for cell-type annotation.

---

## CellGuide Databases

CellGuide is a structured reference database developed by CELLxGENE that **links cell types to their canonical marker genes** across tissues, organs, and species. It aggregates curated annotations from thousands of high-quality single-cell studies, providing standardized marker sets that can be used for gene-set scoring, cluster interpretation, and validation of cell-type identities.

### 1. Per-Gene Database (`cellguide_marker_gene_database_per_gene.json`)
Answers: *Which cell types is this gene a marker for?*

```python
marker_db_per_gene["Gad1"]["Mus musculus"]["brain"]
# → ["inhibitory neuron", "GABAergic neuron", ...]
```

**Use case:** Automated cell type annotation.  
**Ranking:** Ordered by marker strength (first = strongest).  
**Preferred for:** Programmatic lookup in annotation pipelines.

### 2. Per-CellType Database (`cellguide_marker_gene_database_per_celltype.json`)
Answers: *What are the marker genes for this cell type?*

```python
marker_db_per_celltype["Mus musculus"]["Brain"]["inhibitory neuron"]
# → ["Gad1", "Gad2", "Slc32a1", ...]
```

**Use case:** Validation and manual annotation.  
**Ranking:** Genes ordered by marker strength.

### 3. Tissue Options (`tissues_per_organism.json`)
Answers: *What tissues are available for this organism?*

```python
tissues_per_organism["Mus musculus"]
# → ["brain", "liver", "kidney", "heart", ...]
```

**Use case:** Input validation or UI dropdown for selecting tissue context.  
Ensures you only query supported organism/tissue combinations.

### Download Databases from Latch Data
All databases can be loaded directly from **Latch Data** using the `LPath` API.

```python
import json
from latch.ldata.path import LPath
from pathlib import Path

def load_json_lpath(uri):
    lp = LPath(uri)
    local = Path(f"{lp.node_id()}.json")
    lp.download(local, cache=True)
    with open(local) as f:
        return json.load(f)

marker_db_per_gene = load_json_lpath("latch:///cellguide_marker_gene_database_per_gene.json")
marker_db_per_celltype = load_json_lpath("latch:///cellguide_marker_gene_database_per_celltype.json")
tissues_per_organism = load_json_lpath("latch:///tissues_per_organism.json")

# Default database for automated annotation
db_path = Path("cellguide_marker_gene_database_per_gene.json")
```
