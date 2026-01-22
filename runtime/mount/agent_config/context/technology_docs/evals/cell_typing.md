Load `technology_docs/evals/cell_types.json` priors and select matching entry by (organism, tissue, context).

Using priors as ground truth, systematically verify all 6 criteria below: 
1. **Major types present**: All major cell types in final annotations must use names from the priors vocabulary exactly. 
2. **Proportion sanity**: Rare types are not dominant; common types are present at expected frequencies.
3. **Marker enrichment**: Each labeled type has ≥2–3 enriched positive markers (primary or alternate) for that cell type.
4. **Anti-marker check**: Anti-markers are not enriched; if violated, flag as mislabel.
5. **Spatial patterns**: Use `capture_widget_image` to verify patterns match priors 
5. **Unknowns**: Percentage of unknown labels <5%. 
