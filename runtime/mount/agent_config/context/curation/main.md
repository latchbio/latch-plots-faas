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

<pre_analysis_questions>
- Do you have a GSE ID to download from GEO?
- Do you have paper text (abstract, methods, results)?
- Do you have an H5AD or need to construct from raw files?
- What organism (human, mouse, other)?
- Need to map cell types to Cell Ontology?
</pre_analysis_questions>

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
- `latch_sample_id`: Unique sample identifier (required)
- `latch_cell_type_lvl_1`: Cell Ontology term ("name/CL:NNNNNNN")
- `latch_disease`: MONDO term ("name/MONDO:NNNNNNN")
- `latch_tissue`: UBERON term ("name/UBERON:NNNNNNN")
- `latch_organism`: e.g., "Homo sapiens", "Mus musculus"
- `latch_sequencing_platform`: EFO term ("name/EFO:NNNNNNN")
- `latch_subject_id`: Subject/donor identifier (if applicable)
- `latch_condition`: Experimental condition (free text)
- `latch_sample_site`: Collection site (if applicable)
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
- No "unknown" where paper provides information
- Cell type proportions biologically plausible
</self_eval_criteria>
