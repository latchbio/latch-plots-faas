# Cell Type Annotation Guide 

## Critical Decision (MUST ASK)
Ask user:
1. "Are you confident in your current clustering quality?"
2. "Do you prefer fast exploratory annotation or slower detailed annotation?"

**Choose approach:**
- Unsatisfactory clusters OR want fast (~minutes per sample) → **Gene Set Scoring** (`technology_docs/cell_type_annotation/gene_set_scoring.md`)
- Good clusters AND can wait >30 minutes → **Cluster-based** (`technology_docs/cell_type_annotation/cluster_based.md`)

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
