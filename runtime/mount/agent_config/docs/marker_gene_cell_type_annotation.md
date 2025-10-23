# Marker Gene-Based Cell Type Annotation

## Purpose

This workflow provides **automated cell type annotation** using the **CellGuide marker gene database**. After identifying marker genes via differential expression, this method queries a curated database to suggest cell type identities based on marker consensus.

## When to Use

- **After differential gene expression analysis** — automatically run this workflow to suggest cell types for each cluster
- **When you have marker genes** — works with any ranked list of cluster-specific genes
- **For data-driven annotation** — leverages a comprehensive, organism- and tissue-aware marker database

## Database Overview

The CellGuide marker gene database is structured as:

```python
{
    "gene_name": {
        "organism_name": {
            "tissue_name": ["cell_type_1", "cell_type_2", ...],
            "All Tissues": ["cell_type_1", "cell_type_2", ...]
        }
    }
}
```

**Key features:**

- **Organism-specific**: e.g., "Mus musculus", "Homo sapiens"
- **Tissue-specific**: e.g., "Brain", "Liver", "All Tissues"
- **Ranked cell types**: order matters — earlier entries are stronger markers
- **Cell ontology terms**: includes varying granularity (e.g., "neuron" vs "excitatory neuron")

### Tissue Options Database

A companion file `tissues_per_organism.json` provides valid tissue options per organism:

```python
{
    "Mus musculus": ["All Tissues", "Brain", "Liver", "Kidney", ...],
    "Homo sapiens": ["All Tissues", "Brain", "Liver", "Heart", ...]
}
```

This file is used to **dynamically populate** the tissue dropdown based on the selected organism.

---

## Workflow

### Step 1: Load the CellGuide Database and Tissue Options

Download and load both the marker gene database and tissue options from Latch Data:

```python
import json
from latch.ldata.path import LPath
from pathlib import Path

# Download marker gene database
db_lpath = LPath("latch:///cellguide_marker_gene_database.json")
local_db_path = Path(f"{db_lpath.node_id()}.json")
db_lpath.download(local_db_path, cache=True)

# Load marker database into memory
with open(local_db_path, "r") as f:
    marker_db = json.load(f)

# Download tissue options database
tissues_lpath = LPath("latch:///tissues_per_organism.json")
local_tissues_path = Path(f"{tissues_lpath.node_id()}.json")
tissues_lpath.download(local_tissues_path, cache=True)

# Load tissue options into memory
with open(local_tissues_path, "r") as f:
    tissues_per_organism = json.load(f)
```

### Step 2: Determine Organism and Tissue Context

Ask the user to specify organism and tissue type. This context is **critical** for accurate annotation.

The tissue dropdown is **dynamically populated** based on the selected organism using the `tissues_per_organism` database.

```python
from lplots.widgets.select import w_select

# Organism selection
organism = w_select(
    label="Organism",
    options=list(tissues_per_organism.keys()),
    default="Mus musculus" if "Mus musculus" in tissues_per_organism else list(tissues_per_organism.keys())[0]
)

# Get available tissues for selected organism
available_tissues = tissues_per_organism.get(organism.value, ["All Tissues"])

# Tissue selection (dynamically populated based on organism)
tissue = w_select(
    label="Tissue type",
    options=available_tissues,
    default="All Tissues" if "All Tissues" in available_tissues else available_tissues[0]
)
```

### Step 3: Query Database for Marker Genes

For each cluster, extract the **top 10 marker genes** from differential expression results and query the database.

**Tissue selection logic:**

1. For each gene, check available tissues: `marker_db[gene][organism].keys()`
2. If the user-specified tissue exists, use it
3. Otherwise, fall back to `"All Tissues"` if available
4. If neither exists, skip the gene

```python
import scanpy as sc
import pandas as pd
from collections import Counter

def query_marker_database(
    marker_db: dict,
    gene_list: list,
    organism: str,
    tissue: str,
    top_n: int = 10
) -> dict:
    """
    Query the CellGuide database for cell type suggestions.

    Returns:
        dict: {cell_type: {"score": int, "supporting_genes": [genes]}}
    """
    cell_type_votes = {}

    for gene in gene_list[:top_n]:
        if gene not in marker_db:
            continue

        if organism not in marker_db[gene]:
            continue

        # Select tissue
        tissues_available = marker_db[gene][organism].keys()
        if tissue in tissues_available:
            selected_tissue = tissue
        elif "All Tissues" in tissues_available:
            selected_tissue = "All Tissues"
        else:
            continue

        # Get ranked cell types for this gene
        cell_types = marker_db[gene][organism][selected_tissue]

        # Score by rank (higher score = better marker)
        # First position gets highest score, decreasing linearly
        for rank, cell_type in enumerate(cell_types):
            score = len(cell_types) - rank  # Position-based score

            if cell_type not in cell_type_votes:
                cell_type_votes[cell_type] = {
                    "score": 0,
                    "supporting_genes": []
                }

            cell_type_votes[cell_type]["score"] += score
            cell_type_votes[cell_type]["supporting_genes"].append(gene)

    return cell_type_votes

# Example: Query for each cluster
organism_val = organism.value
tissue_val = tissue.value

cluster_annotations = {}

for cluster in adata.obs['leiden'].cat.categories:
    # Get marker genes for this cluster
    marker_df = sc.get.rank_genes_groups_df(adata, group=cluster)
    gene_list = marker_df['names'].tolist()

    # Query database
    suggestions = query_marker_database(
        marker_db,
        gene_list,
        organism_val,
        tissue_val,
        top_n=10
    )

    # Sort by score
    sorted_suggestions = sorted(
        suggestions.items(),
        key=lambda x: x[1]["score"],
        reverse=True
    )

    cluster_annotations[cluster] = sorted_suggestions
```

### Step 4: Display Annotation Results

Present the suggested cell types with confidence scores and supporting evidence.

```python
from lplots.widgets.table import w_table

# Create summary table
annotation_summary = []

for cluster, suggestions in cluster_annotations.items():
    if len(suggestions) > 0:
        # Top suggestion
        top_cell_type = suggestions[0][0]
        top_score = suggestions[0][1]["score"]
        supporting_genes = ", ".join(suggestions[0][1]["supporting_genes"][:5])

        # Alternatives
        alternatives = [s[0] for s in suggestions[1:4]]  # Top 3 alternatives

        annotation_summary.append({
            "Cluster": cluster,
            "Suggested Cell Type": top_cell_type,
            "Confidence Score": top_score,
            "Supporting Markers": supporting_genes,
            "Alternative Suggestions": ", ".join(alternatives) if alternatives else "None"
        })
    else:
        annotation_summary.append({
            "Cluster": cluster,
            "Suggested Cell Type": "Unknown",
            "Confidence Score": 0,
            "Supporting Markers": "No matches",
            "Alternative Suggestions": "None"
        })

annotation_df = pd.DataFrame(annotation_summary)
w_table(source=annotation_df, label="Cell Type Annotation Results")
```

### Step 5: Detailed Evidence per Cluster

Show detailed supporting evidence for each cluster's top suggestions.

```python
# Display detailed evidence for each cluster
for cluster, suggestions in cluster_annotations.items():
    if len(suggestions) == 0:
        continue

    # Create detailed table for top 5 suggestions
    detail_rows = []
    for cell_type, info in suggestions[:5]:
        detail_rows.append({
            "Cell Type": cell_type,
            "Score": info["score"],
            "Supporting Genes": ", ".join(info["supporting_genes"]),
            "Gene Count": len(info["supporting_genes"])
        })

    detail_df = pd.DataFrame(detail_rows)
    w_table(
        source=detail_df,
        label=f"Cluster {cluster}: Top Cell Type Candidates"
    )
```

### Step 6: Visualize Annotations

Add the suggested cell types to the AnnData object and visualize on UMAP and spatial embeddings.

```python
from lplots.widgets.h5 import w_h5

# Add top annotation to adata.obs
adata.obs['predicted_cell_type'] = adata.obs['leiden'].map(
    {cluster: suggestions[0][0] if len(suggestions) > 0 else "Unknown"
     for cluster, suggestions in cluster_annotations.items()}
)

# Add confidence scores
adata.obs['annotation_confidence'] = adata.obs['leiden'].map(
    {cluster: suggestions[0][1]["score"] if len(suggestions) > 0 else 0
     for cluster, suggestions in cluster_annotations.items()}
)

# Visualize with w_h5
viewer = w_h5(ann_data=adata)
```

### Step 7: Validate with Marker Expression

Optionally, validate the suggestions by computing gene set scores for the predicted cell types.

```python
import scanpy as sc

# For each unique predicted cell type, compute marker scores
unique_cell_types = adata.obs['predicted_cell_type'].unique()

for cell_type in unique_cell_types:
    if cell_type == "Unknown":
        continue

    # Collect all genes that voted for this cell type
    marker_genes = []
    for cluster, suggestions in cluster_annotations.items():
        for ct, info in suggestions:
            if ct == cell_type:
                marker_genes.extend(info["supporting_genes"])
                break

    # Remove duplicates and filter to genes present in dataset
    marker_genes = list(set(marker_genes))
    marker_genes = [g for g in marker_genes if g in adata.var_names]

    if len(marker_genes) > 0:
        score_name = f"{cell_type}_validation_score"
        sc.tl.score_genes(adata, marker_genes, score_name=score_name)
```

---

## Handling Edge Cases

### Missing Genes

If a marker gene is not in the database:

- **Skip it** and continue with the remaining markers
- If **no markers** match, label cluster as "Unknown"
- Display a warning in the output

### Organism Mismatch

Always confirm organism with the user before querying. If the database lacks the organism:

- Notify the user
- Fall back to manual annotation workflow

### Tissue Mismatch

Use this priority order:

1. **User-specified tissue** (if available in database)
2. **"All Tissues"** (if available)
3. **Skip gene** (if neither exists)

### Similar Cell Types at Different Granularities

The database may return multiple related terms (e.g., "neuron", "excitatory neuron", "cortical excitatory neuron").

**Recommendation:**

- **Show top 3-5 suggestions** to capture different granularities
- **Prefer more specific terms** when confidence is high
- **Use broader terms** when evidence is weak or conflicting

### Low Confidence Annotations

If the top suggestion has a **low score** (e.g., < 10):

- Flag as **"Low Confidence"**
- Recommend manual validation
- Suggest running marker scoring or literature review

---

## Integration with Workflows

### CosMX / Spatial Transcriptomics

This marker gene annotation workflow should **automatically run** after:

1. Differential gene expression analysis (Step 4 in CosMX workflow)
2. Before manual cell type annotation (Step 5 in CosMX workflow)

The database-driven approach provides:

- **Initial suggestions** to guide manual annotation
- **Validation** of manually assigned cell types
- **Discovery** of unexpected cell populations

### Usage Pattern

```python
# 1. Run differential expression (from CosMX workflow)
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

# 2. Automatically run marker database annotation (this workflow)
# ... (Steps 1-6 above) ...

# 3. Optionally refine with manual marker scoring (from CosMX workflow)
# ... (manual validation if needed) ...
```

---

## Best Practices

1. **Always confirm organism and tissue** before querying
2. **Use top 10 markers** per cluster for consensus (balance coverage and noise)
3. **Display supporting evidence** so users can evaluate suggestions
4. **Show multiple alternatives** to account for ontology granularity
5. **Validate with expression** — compute marker scores to verify predictions
6. **Flag low-confidence** annotations for manual review
7. **Cache the database** — download once per session for efficiency

---

## Summary

This workflow provides **automated, data-driven cell type annotation** by:

- Querying a comprehensive marker gene database
- Aggregating evidence across multiple markers
- Ranking suggestions by consensus and marker strength
- Providing transparent, interpretable results

Use this as a **first-pass annotation** to guide downstream analysis and manual refinement.
