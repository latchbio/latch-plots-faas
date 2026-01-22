# Curation Pipeline

Standardizes molecular data for publication and sharing. Ensures data adheres to controlled vocabularies (Cell Ontology, MONDO, UBERON) with required metadata fields.

<detection>
**When to use:**
- User mentions "curate", "harmonize", "standardize", or "publish"
- Working with external data (paper, collaborator, GSE)
- Has paper text to incorporate into analysis
- Needs metadata cleanup for sharing

**Platform integration:**
Interleave with platform analysis (Xenium, Takara, Vizgen):
- Download: REPLACE (no existing data) or run ALONGSIDE platform loading
- Construct Counts: AFTER download, BEFORE platform QC
- QC & Transform: SUPPLEMENT platform steps with paper context
- Cell typing: AUGMENT with ontology mapping and paper context
- Metadata harmonization: ADD to existing analysis
- Validation: FINAL step before export
</detection>

<inputs>
Only ask for:
1. **GSE ID or uploaded files** (required) - Starting point for data acquisition
2. **Paper text** (strongly encouraged) - Abstract, methods, results for context
</inputs>

<inference>
Infer everything else from data - do not ask upfront:

| Field | Inference Source |
|-------|------------------|
| Organism | Gene ID prefix: ENSG→human, ENSMUSG→mouse |
| File format | Downloaded file extensions (.h5ad, .mtx, .csv) |
| Sample structure | GEO metadata `title`, `source_name_ch1` columns |
| Cell counts | Paper text mentions (e.g., "10,234 cells") |
| Cell types | Clustering + marker genes → Cell Ontology |
| Tissue | GEO `source_name_ch1` or paper methods |
| Disease | GEO characteristics or paper title/methods |
| Sequencing platform | SRA metadata `instrument` field |
| Sex/Age | GEO sample characteristics |

**Only ask user when:**
- Inference is ambiguous (multiple valid interpretations)
- Validation fails and needs clarification
- Critical metadata missing from all sources
</inference>

<plan>
1. Download -> `steps/download.md`
2. Construct Counts -> `steps/construct_counts.md`
3. QC & Filtering -> `steps/qc_filtering.md`
4. Transform & Integration -> `steps/transform_integration.md`
5. Cell Type Annotation -> `steps/cell_typing.md`
6. Metadata Harmonization -> `steps/harmonize_metadata.md`
7. Final Validation & Export -> `steps/validation.md`
</plan>

<interleave_guidance>
**Download:** Run parallel with platform loading, or to retrieve data. Paper text benefits all downstream steps.

**Construct Counts:** Run if parsing raw GEO files. Skip if H5AD exists.

**QC/Transform:** Skip if platform QC comprehensive, or add LLM-guided adaptive thresholds. Ensure `latch_sample_id` exists.

**Cell typing:** Map existing labels to Cell Ontology. Use paper context to refine.

**End of analysis:** Add metadata harmonization + final validation.
</interleave_guidance>

<required_obs_columns>
| Column | Format | Inference Source |
|--------|--------|------------------|
| `latch_sample_id` | string | GEO title/source or user input |
| `latch_cell_type_lvl_1` | "name/CL:NNNNNNN" | Cell typing results |
| `latch_disease` | "name/MONDO:NNNNNNN" | Paper/GEO metadata |
| `latch_tissue` | "name/UBERON:NNNNNNN" | GEO source_name or paper |
| `latch_organism` | "Homo sapiens" / "Mus musculus" | Gene ID prefix |
| `latch_sequencing_platform` | "name/EFO:NNNNNNN" | SRA instrument field |
| `latch_subject_id` | string | GEO characteristics |
| `latch_condition` | free text | Paper/GEO metadata |
| `latch_sample_site` | string | GEO characteristics |
</required_obs_columns>

<required_var_columns>
- var.index: Ensembl gene IDs (ENSG/ENSMUSG format)
- `gene_symbols`: Human-readable gene names (unique)
</required_var_columns>

<self_eval_criteria>
- Required obs columns present with valid ontology terms
- var.index contains valid Ensembl IDs
- gene_symbols exists and is unique
- Counts are raw integers (ingestion) or normalized (after transform)
- No "unknown" where paper/metadata provides information
- Cell type proportions biologically plausible
- Organism matches gene ID prefix (ENSG=human, ENSMUSG=mouse)
</self_eval_criteria>
