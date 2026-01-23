# Step 3 — Cell Type Annotation

<goal>
Assign biologically meaningful cell-type labels to clusters using DGE markers, guided by tissue-specific priors when available, otherwise falling back to CellGuide.
</goal>

<self_eval_criteria>
- Follow `technology_docs/evals/cell_typing.md`
</self_eval_criteria>

<method>

## 1. Check for Tissue-Specific Priors

First, check if `technology_docs/evals/cell_types.json` contains an entry matching the organism + tissue:

- Load and parse `cell_types.json`
- Search `entries` for a match where `match.organism` contains the organism AND `match.tissue` contains the tissue
- If a match is found, use the `priors` from that entry (expected_major_types, expected_subtypes, marker_priors, spatial_morphology_priors)
- If no match, fall back to CellGuide (Section 2)

When priors are available:
- Use `marker_priors` to validate cluster marker genes against known cell type markers
- Use `expected_major_types` and `expected_subtypes` to constrain/prioritize annotations
- Use `anti_markers` to rule out incorrect assignments
- Store the matched prior entry ID in `adata.uns["cell_type_priors_id"]`

## 2. CellGuide Fallback

If no tissue-specific priors exist:

- Load the CellGuide marker database from `latch:///cellguide_marker_gene_database_per_gene.json`
- Use `summarize_clusters()` to match cluster markers to CellGuide cell types

## 3. Assign Labels

- For each cluster, determine the best cell type assignment based on marker overlap
- Write final labels to `adata.obs["cell_type"]`
- Store the full annotation summary in `adata.uns["cell_type_annotation"]`

</method>

<prerequisites>
- Cluster assignments in `adata.obs` (e.g., `clustering_leiden_0.4`)
- DGE results in `adata.uns` (e.g., `rank_genes_groups_clustering_leiden_0.4`)
- Organism and tissue strings (e.g., `"Mus musculus"`, `"brain"`)
</prerequisites>

<library>
- `scanpy`
- `json`
- `xenium_cell_type` helpers (`load_json_lpath`, `summarize_clusters`, etc.)
- Priors: `/opt/latch/plots-faas/runtime/mount/agent_config/context/technology_docs/evals/cell_types.json`
- CellGuide: `latch:///cellguide_marker_gene_database_per_gene.json`
</library>
