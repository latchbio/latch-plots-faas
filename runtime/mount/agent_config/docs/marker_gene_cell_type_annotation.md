# Marker Gene-Based Cell Type Annotation

## Purpose

This workflow provides **automated cell type annotation** using the **CellGuide marker gene database**. After identifying marker genes via differential expression, this method queries a curated database to suggest cell type identities based on marker consensus.

## When to Use

- **After differential gene expression analysis** — automatically run this workflow to suggest cell types for each cluster
- **When you have marker genes** — works with any ranked list of cluster-specific genes
- **For data-driven annotation** — leverages a comprehensive, organism- and tissue-aware marker database

## Database Overview

CellGuide provides **two complementary databases** that enable different types of queries:

### 1. Per-Gene Database (`cellguide_marker_gene_database_per_gene.json`)

**Structure:**

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

**Use case:** "Which cell types is this gene a marker for?"

- **Primary purpose**: Automated cell type annotation from marker genes
- **Workflow**: Given a list of marker genes → find candidate cell types
- **Ranking**: Cell types are ranked by marker strength (first = best)

**Example query:**

```python
marker_db_per_gene["Gad1"]["Mus musculus"]["brain"]
# Returns: ["inhibitory neuron", "GABAergic neuron", ...]
```

### 2. Per-CellType Database (`cellguide_marker_gene_database_per_celltype.json`)

**Structure:**

```python
{
    "organism_name": {
        "tissue_name": {
            "cell_type_name": ["gene_1", "gene_2", ...],
            "All Tissues": {
                "cell_type_name": ["gene_1", "gene_2", ...]
            }
        }
    }
}
```

**Use case:** "What are the marker genes for this cell type?"

- **Primary purpose**: Validation, manual annotation, and marker lookup
- **Workflow**: Given a cell type → find its canonical marker genes
- **Ranking**: Genes are ranked by marker strength (first = best)

**Example query:**

```python
marker_db_per_celltype["Mus musculus"]["Brain"]["inhibitory neuron"]
# Returns: ["Gad1", "Gad2", "Slc32a1", ...]
```

### Key Features (Both Databases)

- **Organism-specific**: e.g., "Mus musculus", "Homo sapiens"
- **Tissue-specific**: e.g., "Brain", "Liver", "All Tissues"
- **Ranked lists**: order matters — earlier entries are stronger markers
- **Cell ontology terms**: includes varying granularity (e.g., "neuron" vs "excitatory neuron")

### When to Use Each Database

| Task                                   | Database to Use  |
| -------------------------------------- | ---------------- |
| Automated annotation from marker genes | **Per-gene**     |
| Validation of predicted cell types     | **Per-celltype** |
| Manual cell type lookup                | **Per-celltype** |
| "What cell types express gene X?"      | **Per-gene**     |
| "What genes mark cell type Y?"         | **Per-celltype** |

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

### Important: Automatic Execution Pattern

**All steps in this workflow should execute automatically** without requiring manual button clicks or other user-triggered actions.

- ✅ **Use widgets for user INPUT**: `w_select`, `w_text_input`, etc. to collect parameters (organism, tissue)
- ❌ **Do NOT use buttons to trigger execution**: avoid patterns like "Click to run annotation" buttons

**Why this matters:**

If essential operations (loading databases, running annotation) are placed behind manual triggers:

1. Later cells will expect the results to exist
2. If the user hasn't clicked the button, those cells will fail
3. The agent may write duplicate auto-running code, leaving dead code in the notebook

**Correct pattern:**

```python
# Collect user input with widgets
organism = w_select(label="Organism", options=[...])
tissue = w_select(label="Tissue", options=[...])

# Then immediately run the annotation (no button needed)
annotations = run_annotation(organism.value, tissue.value)
```

**Incorrect pattern:**

```python
# ❌ Don't do this
organism = w_select(label="Organism", options=[...])
run_button = w_button(label="Run Annotation")

if run_button.clicked:  # This creates dependencies on manual interaction
    annotations = run_annotation(...)
```

---

### Step 1: Load the CellGuide Databases and Tissue Options

Download and load **both CellGuide databases** plus tissue options from Latch Data:

```python
import json
from latch.ldata.path import LPath
from pathlib import Path

# Download per-gene database (for annotation)
db_per_gene_lpath = LPath("latch:///cellguide_marker_gene_database_per_gene.json")
local_db_per_gene_path = Path(f"{db_per_gene_lpath.node_id()}.json")
db_per_gene_lpath.download(local_db_per_gene_path, cache=True)

# Load per-gene database into memory
with open(local_db_per_gene_path, "r") as f:
    marker_db_per_gene = json.load(f)

# Download per-celltype database (for validation/lookup)
db_per_celltype_lpath = LPath("latch:///cellguide_marker_gene_database_per_celltype.json")
local_db_per_celltype_path = Path(f"{db_per_celltype_lpath.node_id()}.json")
db_per_celltype_lpath.download(local_db_per_celltype_path, cache=True)

# Load per-celltype database into memory
with open(local_db_per_celltype_path, "r") as f:
    marker_db_per_celltype = json.load(f)

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

### Step 3: Query Per-Gene Database for Cell Type Suggestions

For each cluster, extract the **top 10 marker genes** from differential expression results and query the **per-gene database** to find candidate cell types.

**Why use per-gene database here:** We have marker genes and want to find which cell types they represent → use `marker_db_per_gene`.

**Tissue selection logic:**

1. For each gene, check available tissues: `marker_db_per_gene[gene][organism].keys()`
2. If the user-specified tissue exists, use it
3. Otherwise, fall back to `"All Tissues"` if available
4. If neither exists, skip the gene

```python
import scanpy as sc
import pandas as pd
from collections import Counter

def query_marker_database(
    marker_db_per_gene: dict,
    gene_list: list,
    organism: str,
    tissue: str,
    top_n: int = 10
) -> dict:
    """
    Query the per-gene CellGuide database for cell type suggestions.

    Args:
        marker_db_per_gene: The per-gene database (gene -> organism -> tissue -> cell types)
        gene_list: List of marker genes to query
        organism: Organism name (e.g., "Mus musculus")
        tissue: Tissue name (e.g., "Brain")
        top_n: Number of top marker genes to use

    Returns:
        dict: {cell_type: {"score": int, "supporting_genes": [genes]}}
    """
    cell_type_votes = {}

    for gene in gene_list[:top_n]:
        if gene not in marker_db_per_gene:
            continue

        if organism not in marker_db_per_gene[gene]:
            continue

        # Select tissue
        tissues_available = marker_db_per_gene[gene][organism].keys()
        if tissue in tissues_available:
            selected_tissue = tissue
        elif "All Tissues" in tissues_available:
            selected_tissue = "All Tissues"
        else:
            continue

        # Get ranked cell types for this gene
        cell_types = marker_db_per_gene[gene][organism][selected_tissue]

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

    # Query per-gene database for cell type suggestions
    suggestions = query_marker_database(
        marker_db_per_gene,  # Use per-gene database
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

### Step 7: Validate with Per-CellType Database

Optionally, validate the suggestions by looking up **canonical marker genes** from the per-celltype database and computing expression scores.

**Why use per-celltype database here:** We have predicted cell types and want to find their expected marker genes → use `marker_db_per_celltype`.

**Two validation approaches:**

#### Approach 1: Validate using genes from annotation (simpler)

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

#### Approach 2: Validate using canonical markers from per-celltype database (recommended)

This approach validates predictions against the **canonical marker gene set** from CellGuide.

```python
import scanpy as sc

# For each unique predicted cell type, lookup canonical markers
unique_cell_types = adata.obs['predicted_cell_type'].unique()

for cell_type in unique_cell_types:
    if cell_type == "Unknown":
        continue

    # Lookup canonical markers from per-celltype database
    canonical_markers = None

    # Try user-specified tissue first
    if organism_val in marker_db_per_celltype:
        if tissue_val in marker_db_per_celltype[organism_val]:
            canonical_markers = marker_db_per_celltype[organism_val][tissue_val].get(cell_type)

        # Fallback to "All Tissues" if not found
        if canonical_markers is None and "All Tissues" in marker_db_per_celltype[organism_val]:
            canonical_markers = marker_db_per_celltype[organism_val]["All Tissues"].get(cell_type)

    if canonical_markers is None:
        continue

    # Filter to genes present in dataset (use top 10 canonical markers)
    canonical_markers_present = [g for g in canonical_markers[:10] if g in adata.var_names]

    if len(canonical_markers_present) > 0:
        score_name = f"{cell_type}_canonical_score"
        sc.tl.score_genes(adata, canonical_markers_present, score_name=score_name)
```

**Comparing validation scores:**

```python
# Compare canonical marker scores across clusters
score_cols = [col for col in adata.obs.columns if col.endswith("_canonical_score")]

if score_cols:
    validation_summary = adata.obs[score_cols + ["leiden", "predicted_cell_type"]].groupby("leiden").mean()
    w_table(source=validation_summary, label="Canonical Marker Expression Scores by Cluster")
```

**Interpretation:**

- **High canonical scores** in predicted cluster → Strong validation
- **Low canonical scores** → Prediction may be incorrect or cell type is ambiguous
- **Mismatched scores** (high score in different cluster) → May indicate annotation error

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

## Database Usage Patterns

Understanding **when to use each database** is critical for efficient and accurate annotation workflows.

### When to Use Per-Gene Database

**Use `marker_db_per_gene`** when you have genes and need to find cell types:

| Scenario                    | Example                                                                                 |
| --------------------------- | --------------------------------------------------------------------------------------- |
| **Automated annotation**    | You have marker genes from differential expression → query to find cell type candidates |
| **Gene-to-celltype lookup** | User asks: "What cell types express Gad1?"                                              |
| **Cluster interpretation**  | Given top genes in a cluster → identify the cell type                                   |

**Code pattern:**

```python
# Given: list of marker genes
marker_genes = ["Gad1", "Gad2", "Slc32a1"]

# Find: which cell types these genes mark
for gene in marker_genes:
    cell_types = marker_db_per_gene[gene][organism][tissue]
    print(f"{gene} marks: {cell_types[:3]}")  # Top 3 cell types
```

### When to Use Per-CellType Database

**Use `marker_db_per_celltype`** when you have cell types and need to find genes:

| Scenario                    | Example                                                                           |
| --------------------------- | --------------------------------------------------------------------------------- |
| **Validation**              | You predicted "inhibitory neuron" → lookup canonical markers to verify expression |
| **Manual annotation**       | User wants to manually check markers for a specific cell type                     |
| **Celltype-to-gene lookup** | User asks: "What are the markers for astrocytes?"                                 |
| **Reference gene sets**     | Building gene sets for scoring or visualization                                   |

**Code pattern:**

```python
# Given: predicted cell type
predicted_cell_type = "inhibitory neuron"

# Find: canonical marker genes for this cell type
marker_genes = marker_db_per_celltype[organism][tissue][predicted_cell_type]
print(f"Top markers for {predicted_cell_type}: {marker_genes[:10]}")  # Top 10 markers
```

### Complete Workflow Example: Using Both Databases

**Step 1:** Annotation (per-gene) → **Step 2:** Validation (per-celltype)

```python
# Step 1: Annotation using per-gene database
# We have marker genes from differential expression
top_cluster_markers = ["Gad1", "Gad2", "Slc32a1", "Pvalb"]

# Query per-gene database to find cell types
cell_type_suggestions = {}
for gene in top_cluster_markers:
    if gene in marker_db_per_gene and organism in marker_db_per_gene[gene]:
        cell_types = marker_db_per_gene[gene][organism].get(tissue, [])
        for ct in cell_types[:5]:  # Top 5 cell types per gene
            cell_type_suggestions[ct] = cell_type_suggestions.get(ct, 0) + 1

# Top suggestion by vote count
predicted_cell_type = max(cell_type_suggestions, key=cell_type_suggestions.get)
print(f"Predicted: {predicted_cell_type}")

# Step 2: Validation using per-celltype database
# Lookup canonical markers for the predicted cell type
canonical_markers = marker_db_per_celltype[organism][tissue].get(predicted_cell_type, [])
print(f"Expected markers: {canonical_markers[:10]}")

# Compute expression score for validation
canonical_markers_in_data = [g for g in canonical_markers[:10] if g in adata.var_names]
sc.tl.score_genes(adata, canonical_markers_in_data, score_name=f"{predicted_cell_type}_validation")
```

### Decision Tree: Which Database?

```
Do you have...
├─ Genes and need cell types? → Use per-gene database
│   ├─ Marker genes from DE? → Automated annotation workflow
│   └─ Single gene lookup? → Direct query: marker_db_per_gene[gene][org][tissue]
│
└─ Cell types and need genes? → Use per-celltype database
    ├─ Predicted cell type to validate? → Canonical marker validation
    └─ Manual cell type to explore? → Direct query: marker_db_per_celltype[org][tissue][cell_type]
```

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
2. **Load both databases** at the start — you'll need both for annotation and validation
3. **Use per-gene for annotation** — when going from markers to cell types
4. **Use per-celltype for validation** — when verifying predictions with canonical markers
5. **Use top 10 markers** per cluster for consensus (balance coverage and noise)
6. **Display supporting evidence** so users can evaluate suggestions
7. **Show multiple alternatives** to account for ontology granularity
8. **Validate with expression** — compute canonical marker scores to verify predictions
9. **Flag low-confidence** annotations for manual review
10. **Cache the databases** — download once per session for efficiency

---

## Summary

This workflow provides **automated, data-driven cell type annotation** using **two complementary databases**:

**Per-Gene Database (`marker_db_per_gene`):**

- Maps genes → cell types
- Used for automated annotation from marker genes
- Primary workflow: differential expression → cell type suggestions

**Per-CellType Database (`marker_db_per_celltype`):**

- Maps cell types → genes
- Used for validation and marker lookup
- Primary workflow: predicted cell type → canonical marker validation

**Complete annotation workflow:**

1. Load both databases and tissue options
2. Use **per-gene database** to annotate clusters from marker genes
3. Use **per-celltype database** to validate predictions with canonical markers
4. Visualize and interpret results

The dual-database approach provides:

- **Automated suggestions** from observed markers
- **Validation** against canonical marker gene sets
- **Transparent evidence** for all predictions
- **Flexibility** for both automated and manual workflows

Use this as a **first-pass annotation** to guide downstream analysis and manual refinement.
