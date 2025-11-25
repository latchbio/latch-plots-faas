# Cell Type Annotation Guide

## Supported Methods

**Per-cluster Gene Set Scoring**  
- Scores canonical marker sets per cluster and assigns one label per group.  
Best when clusters are well separated and biologically meaningful.  
- See `technology_docs/atlasxomics/cell_type_annotation/per_cluster_gene_set_scoring.md`

**Per-cell Gene Set Scoring**  
- Rapid, exploratory annotation at single-cell resolution.  
- Useful for capturing heterogeneity or when clustering is preliminary.  
- See `technology_docs/cell_type_annotation/per_cell_gene_set_scoring.md`

**Unbiased Cluster-based Annotation (DE-driven)**  
- Performs differential expression between clusters to identify defining markers, then matches them to known cell types using literature or databases.  
- Helpful when expected cell types are unclear or curated marker sets are incomplete.
- See `technology_docs/atlasxomics/cell_type_annotation/cluster_based_de.md`

## When to Use Each Method

- Clusters are reliable and you want concise, interpretable labels (one cell type label per cluster) → **Per-cluster scoring**
- You want to capture heterogeneity within clusters or need rapid exploration → **Per-cell scoring**
- Tissue composition is uncertain or marker sets are limited. Unbiased marker discovery by comparing between clusters (slow, 30+ minutes for large datasets) → **DE-driven**

## Input: `*_sm_ge.h5ad` (gene activity scores)

## CellGuide Databases

Curated reference from CELLxGENE linking cell types to canonical marker genes across tissues/species. Aggregates annotations from thousands of single-cell studies with standardized marker sets.

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

# Load databases
marker_db_per_gene = load_json_lpath("latch:///cellguide_marker_gene_database_per_gene.json")
marker_db_per_celltype = load_json_lpath("latch:///cellguide_marker_gene_database_per_celltype.json")
tissues_per_organism = load_json_lpath("latch:///tissues_per_organism.json")
```

#### Usage

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
