# Curation Pipeline

Standardizes molecular data for publication and sharing. Ensures data adheres to controlled vocabularies (Cell Ontology, MONDO, UBERON) with required metadata fields.

<detection>
**When to use:**
- User mentions "curate", "harmonize", "standardize", or "publish"
- Working with external data (paper, collaborator, GSE)
- Has paper text to incorporate into analysis
- Needs metadata cleanup for sharing
</detection>

<inputs>
Only ask for:
1. **GSE ID or uploaded files** (required)
2. **Paper text** (strongly encouraged) - via widget, not chat
</inputs>

<plan>
1. Download → `steps/download.md`
2. Construct Counts → `steps/construct_counts.md`
3. **Merge with technology-specific plan** (see below)
4. Harmonize Metadata → `steps/harmonize_metadata.md` (end of analysis)
</plan>

<platform_merge>
After download + construct_counts, detect platform from `tmp/{accession}/study_metadata.txt`:

| Indicator in metadata | Technology docs |
|-----------------------|-----------------|
| "Xenium" or GPL33896 | `technology_docs/xenium/main.md` |
| "Visium" | `technology_docs/visium.md` |
| "Vizgen" or "MERFISH" | `technology_docs/vizgen/main.md` |
| "Takara" or "ICELL8" | `technology_docs/takara/main.md` |
| "AtlasXomics" | `technology_docs/atlasxomics/main.md` |

**Follow the technology-specific plan for QC, normalization, clustering, cell typing.**

At each step, supplement with curation context:
- `study_metadata.txt` → sample info, organism, tissue
- `paper_text.txt` → QC thresholds, expected cell counts, cell types mentioned
</platform_merge>

<inference>
Infer from data - do not ask upfront:

| Field | Source |
|-------|--------|
| Organism | Gene ID prefix (ENSG→human, ENSMUSG→mouse) |
| Platform | `study_metadata.txt` platform_id, label_ch1 |
| Tissue | GEO source_name_ch1 or paper methods |
| Disease | GEO characteristics or paper |
| Cell types | Analysis results → Cell Ontology |

Only ask when inference is ambiguous or validation fails.
</inference>

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
