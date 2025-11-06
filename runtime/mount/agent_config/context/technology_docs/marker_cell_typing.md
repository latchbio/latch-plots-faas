# Marker-based Cell Type Annotation

## Overview

This guide provides a step-by-step pipeline for annotating cell types in spatial biology datasets using marker gene analysis and the CellGuide database.

There are two general strategies depending on data quality and clustering resolution:

- **Approach 1**: Gene set scoring — uses predefined marker gene sets to assign cell types directly; _does not_ require well-separated clusters.

- **Approach 2**: Cluster-based annotation — identifies marker genes per cluster and matches them to known cell types; requires clear, biologically meaningful clusters.

**ALWAYS confirm with the user which approach to use before proceeding.** Do _not_ automatically default to cluster-based annotation just because "cluster" exists in `adata.obs`.

--

## CellGuide Database Setup

### Database Overview

CellGuide provides **two complementary databases** that enable different types of queries:

#### 1. Per-Gene Database (`cellguide_marker_gene_database_per_gene.json`)

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

#### 2. Per-CellType Database (`cellguide_marker_gene_database_per_celltype.json`)

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

#### 3. Tissue Options Database (`tissues_per_organism.json`)

**Structure:**

```python
{
    "organism_name": ["tissue_1", "tissue_2", "tissue_3", ...]
}
```

**Use case:** "What tissues are available for this organism?"

- **Primary purpose**: Provide valid tissue options for a given organism when querying CellGuide databases
- **Workflow**: Given an organism → get list of available tissues for that organism
- **Use in UI**: Populate dropdowns or validate user-selected tissue names

**Example query:**

```python
tissues_per_organism["Mus musculus"]
# Returns: ["brain", "liver", "kidney", "heart", ...]
```

#### Key Features (All Databases)

- **Organism-specific**: e.g., `"Mus musculus"`, `"Homo sapiens"`
- **Tissue-specific**: e.g., `"Brain"`, `"Liver"`, `"All Tissues"` (fallback)
- **Ranked lists**: order matters — earlier entries are stronger markers
- **Cell ontology terms**: includes varying granularity (e.g., `"neuron"` vs `"excitatory neuron"`)

**Note**: This guide uses the **per-gene database** for automated annotation. The per-celltype database is useful for validation and manual annotation workflows. The **tissue options database** helps determine available tissues for a given organism, which is useful for UI validation and user guidance.

### Downloading from Latch Data

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

# For annotation workflow, use the per-gene database path
db_path = local_db_per_gene_path
```

---
# **Approach 1**: Gene set scoring

This approach assigns cell types by scoring predefined marker gene sets across cells, without requiring prior clustering. For multi-sample datasets, scoring should be performed per condition or sample to capture biological variation.

---
# **Approach 2**: Cluster-based Annotation

The workflow identifies marker genes per cluster, finds consensus markers across conditions, and matches them against known cell type markers in CellGuide.

## Prerequisites

- **AnnData object** (`adata`) with:
  - `adata.obs["cluster"]` - cluster assignments for each cell (required)
  - Sample/condition identifier column (optional, varies by dataset):
    - Common names: `"sample"`, `"condition"`, `"batch"`
    - **Must introspect `adata.obs.columns` to determine available columns**
  - Gene expression data in `adata.X` or `adata.layers`
- **CellGuide database** (`cellguide_marker_gene_database_per_gene.json`) - JSON file containing marker genes for known cell types
  - May need to download from Latch Data (see Database Setup below)
- **Organism and tissue context** - Required for accurate CellGuide lookups (e.g., `"Mus musculus"`, `"brain"`)

**Important**: Different `adata` objects use different column names. Always introspect the structure before making assumptions.

## Step-by-step Guide

### Step 0: Introspect AnnData Structure

**Objective**: Identify available columns and determine if the dataset has multiple samples/conditions.

**Process**:

- **Always introspect `adata.obs.columns`** to see what metadata columns are available
- Check for cluster column (required): typically `"cluster"` but may be `"leiden"`, `"louvain"`, etc.
- Check for sample/condition columns (optional): common names include:

  - `"sample"`, `"condition"`, `"batch"`, `"sample_id"`

- Determine if multiple samples exist: if a sample column exists and has >1 unique value
- If no sample column exists or only one unique value, treat as single-sample dataset

**Key check**:

- If sample column exists and has multiple unique values → use per-sample analysis
- Otherwise → use single-sample analysis

### Step 1: Perform Differential Expression Analysis

**Objective**: Identify marker genes for each cluster using statistical testing.

**Process**:

- If `adata` contains **multiple samples/conditions**, first split the `adata` into smaller **per-sample** `adata`s. Then, perform a **per-sample 1-versus-all Wilcoxon rank-sum test** using `scanpy.tl.rank_genes_groups()`.
- This ensures robust marker identification across biological replicates.
- Combine results into a single dataframe with ranked marker gene scores across all samples.

**Key parameters**:

- Use `method='wilcoxon'` for non-parametric testing
- Set `n_genes=300` to get top genes. Filter if needed later.
- Group by `"cluster"` column in `adata.obs`
- If multiple samples exist, include the sample identifier in the grouping

**Output**: A pandas DataFrame with columns:

- `"names"` - gene names
- `"group"` - cluster ID
- `"pvals_adj"` - adjusted p-values
- `"scores"` - test statistic scores
- Sample identifier column (name matches detected column from `adata.obs`, only if multiple samples exist)

### Step 2: Find Consensus Marker Genes

**Objective**: Identify marker genes that are consistently expressed across multiple conditions/samples.

**Process**:

- For each cluster, filter genes with `pvals_adj < 0.05`
- Count how many conditions/samples each gene appears as a marker in
- Identify:
  - **`genes_all`**: Genes that appear as markers in **all** conditions/samples (100% consensus)
  - **`genes_70`**: Genes that appear as markers in **≥70%** of conditions/samples (high consensus)

**Rationale**: Consensus markers are more reliable for cell type annotation than condition-specific markers.

### Step 3: Lookup Cell Types in CellGuide Database

**Objective**: Match identified marker genes to known cell types in the CellGuide database.

**Process**:

- For each cluster's consensus markers, query the CellGuide database
- The database is organized by: `gene → organism → tissue → [list of cell types]`
- Lookup priority:
  1. Tissue-specific cell types (if tissue is known)
  2. "All Tissues" entries (if tissue-specific not found)
  3. Skip if no match for organism

**Required inputs**:

- `organism` - e.g., `"Mus musculus"`, `"Homo sapiens"`
- `tissue` - e.g., `"brain"`, `"liver"`, `"All Tissues"` (fallback)
- `db_path` - Path to CellGuide JSON database file

### Step 4: Summarize Cell Type Annotations

**Objective**: Create a comprehensive summary for each cluster with the most likely cell type(s).

**Process**:

- Count how many marker genes support each cell type
- Identify the cell type(s) with the highest marker gene count
- Collect the specific marker genes that support the top cell type(s)
- Generate a summary dictionary for each cluster

**Output structure**:

```python
summary[cluster] = {
    "genes_all": genes_all,                    # Genes in all conditions
    "genes_70": genes_70,                      # Genes in ≥70% conditions
    "most_common_cell_type": top_cell_types,   # List of top cell type(s)
    "cell_type_counts": cell_type_counts,      # Dict: cell_type → marker count
    "markers_for_most_common_cell_type": top_markers,  # Supporting marker genes for top cell type(s)
}
```

**Note**: When multiple cell types have the same highest marker count, `most_common_cell_type` will be a list containing all tied cell types. See "Resolving Tied Cell Types" below for refinement strategies.

### Step 5: Annotate AnnData Object

**Objective**: Add cell type annotations to `adata.obs` as a new column for downstream analysis and visualization.

**Process**:

- Map each cluster to its most common cell type from the summary
- Handle cases where multiple cell types are tied (join with separator)
- Handle clusters with no annotation (use "Unknown")
- Add a new column `adata.obs["cell_type"]` containing the annotations

**Key considerations**:

- **Multiple cell types**: If a cluster has multiple top cell types, join them with a separator (e.g., `", "`)
- **No annotation**: If `most_common_cell_type` is empty, assign `"Unknown"`
- **Cluster matching**: Ensure cluster IDs match between `summary` keys and `adata.obs["cluster"]` values (convert to strings if needed)

**Output**:

- New column `adata.obs["cell_type"]` with cell type annotations for each cell
- Each cell inherits the annotation of its cluster

## Implementation Guide

### Step 1: Differential Expression Analysis

**Single sample dataset**:

```python
import scanpy as sc
import pandas as pd

PVAL_MAX, LFC_MIN, TOP_N = 0.05, 0.5, 300

sc.tl.rank_genes_groups(
    adata,
    groupby="cluster",
    method="wilcoxon",
    n_genes=None,
    use_raw=False,
)

df = sc.get.rank_genes_groups_df(adata, group=None)

df = df[
    (df.pvals_adj < PVAL_MAX) &
    (df.logfoldchanges > LFC_MIN) &
    (df.scores > SCORE_MIN)
]

df = df.sort_values(["group", "logfoldchanges", "scores"], ascending=[True, False, False])

ranked_genes_df = df.groupby("group", as_index=False).head(TOP_N).reset_index(drop=True)
```

**Multiple samples/conditions dataset**:

```python
import scanpy as sc
import pandas as pd

PVAL_MAX, LFC_MIN, TOP_N = 0.05, 0.5, 300

# Perform per-sample marker identification
all_results = []
for sample in adata.obs["sample_id"].unique():
    adata_subset = adata[adata.obs["sample_id"] == sample].copy()

    sc.tl.rank_genes_groups(
        adata_subset,
        groupby="cluster",
        method="wilcoxon",
        n_genes=None,  # all genes
        use_raw=False,
    )

    df = sc.get_rank_genes_groups_df(adata_subset, group=None)
    df["sample_id"] = sample

    # Filter by significance and effect size
    df = df[(df.pvals_adj < PVAL_MAX) & (df.logfoldchanges > LFC_MIN)]

    # Sort by biological and statistical strength
    df = df.sort_values(["group", "logfoldchanges", "scores"], ascending=[True, False, False])

    # Keep top N per cluster
    df = df.groupby("group", as_index=False).head(TOP_N).reset_index(drop=True)
    all_results.append(df)

# Combine all results
ranked_genes_df = pd.concat(all_results, ignore_index=True)
```

### Step 2 & 3: Find Consensus Markers and Lookup Cell Types

```python
import pandas as pd
import json
from collections import Counter
from pathlib import Path

def find_overlap_genes(cluster, ranked_genes_df):
    """
    Find consensus marker genes for a cluster across conditions.
    
    Args:
        cluster: Cluster ID (int or str)
        ranked_genes_df: DataFrame with columns: "names", "group", "pvals_adj", "sample_id"
    
    Returns:
        tuple: (genes_all, genes_70) - lists of gene names
            - genes_all: markers in all conditions
            - genes_70: markers in ≥70% of conditions
    """
    group = ranked_genes_df[ranked_genes_df["group"] == cluster]
    
    # Filter by significance
    group_filter = group[group["pvals_adj"] < 0.05]
    
    if group_filter.empty:
        return [], []
    
    # Count conditions per gene
    overlap = (
        group_filter.groupby("names")["sample_id"]
        .nunique()
        .reset_index(name="n_conditions")
    )
    
    n_total = group_filter["sample_id"].nunique()
    
    # Find consensus genes
    genes_all = overlap.loc[overlap.n_conditions == n_total, "names"]
    genes_70 = overlap.loc[overlap.n_conditions >= 0.7 * n_total, "names"]
    
    return list(genes_all), list(genes_70)


def lookup_cellguide_celltypes(marker_genes, organism, tissue, db_path):
    """
    Lookup cell types in CellGuide database for given marker genes.
    
    Args:
        marker_genes: List of gene names (str)
        organism: Organism name, e.g., "Mus musculus", "Homo sapiens"
        tissue: Tissue name, e.g., "brain", "liver", or "All Tissues"
        db_path: Path to CellGuide JSON database file
    
    Returns:
        dict: {gene_name: [list of cell types]} for genes found in database
    """
    # Load database (handle both local Path and string paths)
    if isinstance(db_path, (str, Path)):
        db_path = Path(db_path)
    
    with open(db_path, "r") as f:
        db = json.load(f)
    
    results = {}
    for gene in marker_genes:
        gene = gene.strip()
        
        # Navigate database structure: gene → organism → tissue
        org_entry = db.get(gene, {}).get(organism)
        if not org_entry:
            continue
        
        # Try tissue-specific first, then fallback to "All Tissues"
        tissue_entry = (
            org_entry.get(tissue)
            or org_entry.get("All Tissues")
            or []
        )
        
        if tissue_entry:
            results[gene] = tissue_entry
    
    return results


def summarize_clusters(
    adata,
    ranked_genes_df,
    db_path,
    organism="Mus musculus",
    tissue="brain"
):
    """
    Annotate clusters with cell types using marker genes and CellGuide database.
    
    Args:
        adata: AnnData object with cluster assignments in adata.obs["cluster"]
        ranked_genes_df: DataFrame from differential expression analysis
        db_path: Path to CellGuide JSON database
        organism: Organism name for CellGuide lookup
        tissue: Tissue name for CellGuide lookup
    
    Returns:
        tuple: (summary, celltype_markers_dict)
            - summary: dict mapping cluster → annotation details
            - celltype_markers_dict: dict mapping cell_type → list of marker genes
    """
    summary = {}
    celltype_markers_dict = {}
    
    for cluster in sorted(adata.obs["cluster"].unique()):
        # Find consensus markers (Step 2)
        genes_all, genes_70 = find_overlap_genes(cluster, ranked_genes_df)
        
        if not genes_all:
            # No consensus markers found
            summary[cluster] = {
                "genes_all": [],
                "genes_70": list(genes_70),
                "most_common_cell_type": [],
                "cell_type_counts": {},
                "markers_for_most_common_cell_type": [],
            }
            continue
        
        # Lookup cell types in CellGuide (Step 3)
        lookup = lookup_cellguide_celltypes(genes_all, organism, tissue, db_path)
        
        if not lookup:
            # No matches in CellGuide
            summary[cluster] = {
                "genes_all": list(genes_all),
                "genes_70": list(genes_70),
                "most_common_cell_type": [],
                "cell_type_counts": {},
                "markers_for_most_common_cell_type": [],
            }
            continue
        
        # Count genes per cell type
        counter = Counter(
            ct for cts in lookup.values() for ct in cts
        )
        cell_type_counts = dict(
            sorted(counter.items(), key=lambda x: (-x[1], x[0]))
        )
        
        # Identify top cell type(s) - those with highest marker count
        if counter:
            max_count = max(counter.values())
            top_cell_types = [
                ct for ct, n in counter.items() if n == max_count
            ]
        else:
            top_cell_types = []
        
        # Collect marker genes supporting top cell type(s)
        top_markers = [
            g for g, cts in lookup.items()
            if any(ct in top_cell_types for ct in cts)
        ]
        
        # Create summary for this cluster
        summary[cluster] = {
            "genes_all": list(genes_all),
            "genes_70": list(genes_70),
            "most_common_cell_type": top_cell_types,
            "cell_type_counts": cell_type_counts,
            "markers_for_most_common_cell_type": top_markers,
        }
        
        # Update global dict (celltype → all markers across clusters)
        for ct in top_cell_types:
            celltype_markers_dict.setdefault(ct, set()).update(top_markers)
    
    # Convert sets to lists for JSON serialization
    celltype_markers_dict = {
        ct: sorted(list(m)) for ct, m in celltype_markers_dict.items()
    }
    
    return summary, celltype_markers_dict
```

### Annotate AnnData Object

```python
def annotate_adata(adata, cluster_summary, cluster_col="cluster", cell_type_col="cell_type"):
    """
    Add cell type annotations to adata.obs based on cluster summary.
    
    Args:
        adata: AnnData object with cluster assignments
        cluster_summary: Dictionary from summarize_clusters() mapping cluster → annotation details
        cluster_col: Name of cluster column in adata.obs (default: "cluster")
        cell_type_col: Name for new cell type column (default: "cell_type")
    
    Returns:
        AnnData object with new cell_type column added to adata.obs
    """
    # Create mapping from cluster to cell type string
    cluster_to_celltype = {}
    
    for cluster, details in cluster_summary.items():
        cell_types = details.get("most_common_cell_type", [])
        
        if cell_types:
            # Join multiple cell types with separator
            cell_type_str = ", ".join(cell_types)
        else:
            # No annotation found
            cell_type_str = "Unknown"
        
        # Convert cluster to string to match adata.obs format
        cluster_to_celltype[str(cluster)] = cell_type_str
    
    # Map cluster IDs to cell types
    adata.obs[cell_type_col] = adata.obs[cluster_col].astype(str).map(
        cluster_to_celltype
    ).fillna("Unknown")
    
    # Report annotation statistics
    n_annotated = (adata.obs[cell_type_col] != "Unknown").sum()
    n_total = len(adata.obs)
    n_unique_types = adata.obs[cell_type_col].nunique()
    
    print(f"Annotated {n_annotated}/{n_total} cells ({100*n_annotated/n_total:.1f}%)")
    print(f"Found {n_unique_types} unique cell types")
    print(f"Cell type distribution:\n{adata.obs[cell_type_col].value_counts()}")
    
    return adata
```

## Important Considerations

### When to Use This Method

- **Best for**: Well-characterized cell types with established marker genes
- **Requires**: Pre-clustered data and known organism/tissue context
- **Limitations**:
  - May not identify novel or rare cell types
  - Depends on quality and completeness of CellGuide database
  - Tissue-specific annotations require accurate tissue information

### Error Handling

- **Missing cluster column**: Check `"cluster"` in `adata.obs.columns` before proceeding
- **Empty marker lists**: If `genes_all` is empty, consider using `genes_70` or relaxing p-value threshold
- **No CellGuide matches**: May indicate:
  - Incorrect organism/tissue specification
  - Novel cell type not in database
  - Poor marker gene quality
- **Multiple top cell types**: When ties occur, all top cell types are returned (may indicate similar cell lineages)
