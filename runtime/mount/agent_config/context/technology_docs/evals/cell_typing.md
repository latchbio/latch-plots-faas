Load `technology_docs/cell_types.json` priors and select matching entry by (organism, tissue, context).

Using priors as ground truth, verify:  
1. **Major types present**: All expected major types appear in annotations.
2. **Proportion sanity**: Rare types are not dominant; common types are present at expected frequencies.
3. **Marker enrichment**: Each labeled type has ≥2–3 enriched positive markers from the panel.
4. **Anti-marker check**: Anti-markers are not enriched; if violated, flag as mislabel.
5. **Spatial patterns**: Use `capture_widget_image` to verify patterns match priors (e.g., follicle layers, vascular tracks, punctate infiltrates).
6. **Cross-sample consistency**: Major-type calls are consistent across samples unless explained by condition/timepoint.
