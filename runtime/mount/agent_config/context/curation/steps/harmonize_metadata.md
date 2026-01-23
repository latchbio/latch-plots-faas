<goal>
Add required latch_* columns with ontology terms, inferred from analysis results and context files.
</goal>

<method>
Read `tmp/{accession}/study_metadata.txt` and `paper_text.txt`. Map inferred values to ontology terms and populate required obs columns. Run at end of analysis after cell typing is complete.
</method>

<required_columns>
| Column | Format | Source |
|--------|--------|--------|
| `latch_sample_id` | string | GEO title or user input |
| `latch_cell_type_lvl_1` | "name/CL:NNNNNNN" | Cell typing results |
| `latch_disease` | "name/MONDO:NNNNNNN" | Paper/metadata |
| `latch_tissue` | "name/UBERON:NNNNNNN" | GEO source_name or paper |
| `latch_organism` | "Homo sapiens" / "Mus musculus" | Gene ID prefix |
| `latch_sequencing_platform` | "name/EFO:NNNNNNN" | Platform from metadata |
</required_columns>
