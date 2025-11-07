# Vizgen MERFISH Transcriptomics Workflow + Evals

This document summarises the MERFISH imaging-based transcriptomics evaluations under:

- `vizgen/cell_typing`
- `vizgen/differential-expression`
- `vizgen/tissue_composition`

Each evaluation lists its biological motivation, analytic objective, inputs, expected output schema, ground-truth reference, tolerance window, and behaviour under test.

---

## Cell Type Marker Recovery (`vizgen/cell_typing`)

**Shared method**

- **What & Why**: All marker-recovery evals take clustered MERFISH expression, run one-vs-all t-tests (logFC > 0, adj p ≤ 0.05), uppercase/deduplicate genes, and confirm that lineage-specific clusters surface the expected canonical markers—evidence that MERFISH captures true biological programs.
- **Goal**: Output `{"markers_predicted": [...]}` containing the top-k markers for the lineage of interest.
- **Input**: `latch://38438.account/Mouse_Aging_subsetted.h5ad`.
- **Expected Output**: Ranked gene list (k defined per eval).
- **Ground Truth & Tolerance**: Lineage-specific marker panels with Jaccard thresholds that correspond to ≥ 80 % overlap (e.g., ≥ 0.50 for 8-gene sets, ≥ 0.29 for 2-gene sets).
- **Core Test**: Demonstrates that unsupervised clustering plus differential testing recovers the biologically appropriate transcriptomic signatures.

**Eval summaries**

- [Astrocyte reactive markers](cell_typing/mouse_brain_merfish_astrocyte_onevsall_marker_jaccard_v9_ttest_symbols.json) — captures `{GFAP, VIM, C4B, C3, SERPINA3N, CXCL10, IL18, HIF3A}`; top 13; Jaccard ≥ 0.50.
- [Endothelial interferon marker](cell_typing/mouse_brain_merfish_endothelial_onevsall_marker_jaccard_v9_ttest_symbols.json) — requires `XDH`; top 6; Jaccard ≥ 0.17.
- [Inhibitory interneuron panel](cell_typing/mouse_brain_merfish_inhibitory_onevsall_marker_jaccard_v9_ttest_symbols.json) — expects `{SST, PVALB, VIP, LAMP5}`; top 9; Jaccard ≥ 0.44.
- [Macrophage injury markers](cell_typing/mouse_brain_merfish_macrophage_onevsall_marker_jaccard_v9_ttest_symbols.json) — `{LYZ2, EMR1, ITGAM}`; top 8; Jaccard ≥ 0.38.
- [Microglial activation set](cell_typing/mouse_brain_merfish_microglia_onevsall_marker_jaccard_v9_ttest_symbols.json) — `{B2M, TREM2, CCL2, APOE, AXL, ITGAX, CD9, C1QA, C1QC, LYZ2, CTSS}`; top 16; Jaccard ≥ 0.50.
- [Medium spiny neuron identity](cell_typing/mouse_brain_merfish_msn_onevsall_marker_jaccard_v9_ttest_symbols.json) — `{DRD1, DRD2}`; top 7; Jaccard ≥ 0.29.
- [Pan-neuronal lineage recall](cell_typing/mouse_brain_merfish_neuron_onevsall_marker_jaccard_v9_ttest_symbols.json) — `{SST, PVALB, VIP, LAMP5, CALB2, DRD1, DRD2, LHX6, CHAT}`; top 14; Jaccard ≥ 0.53.
- [Oligodendrocyte inflammatory state](cell_typing/mouse_brain_merfish_oligodendrocyte_onevsall_marker_jaccard_v9_ttest_symbols.json) — `{C4B, IL33, IL18}`; top 8; Jaccard ≥ 0.38.
- [OPC lineage fidelity](cell_typing/mouse_brain_merfish_opc_onevsall_marker_jaccard_v9_ttest_symbols.json) — `{PDGFRA, CSPG4}`; top 7; Jaccard ≥ 0.29.
- [Pericyte interferon response](cell_typing/mouse_brain_merfish_pericyte_onevsall_marker_jaccard_v9_ttest_symbols.json) — `{IFIT3, IFITM3}`; top 7; Jaccard ≥ 0.29.
- [T cell infiltration](cell_typing/mouse_brain_merfish_tcell_onevsall_marker_jaccard_v9_ttest_symbols.json) — `{CD3E, TRAC, IL7R}`; top 8; Jaccard ≥ 0.38.

---

## Ageing Differential Expression (`vizgen/differential-expression`)

Each lineage ships two evals:

1. `*_onevsone_*` — rediscover cell types via the mandated hierarchy (Neuron vs Non-neuronal → sub-lineages) before differential testing.  
2. `*_ttest_*` — reuse provided labels and run the same OLD vs JUVENILE comparison.

**Shared method**

- **What & Why**: Measure ageing-induced transcriptional programs by contrasting 90 wk (OLD) vs 4 wk (JUVENILE) cells within each lineage, capturing inflammatory, interferon, or degeneration-associated signals.
- **Goal**: Emit `{"markers_predicted": [...]}` comprising the top |ground_truth| + 5 OLD-upregulated genes.
- **Inputs**: Hierarchical evals use `Mouse_Aging_subsetted.h5ad`; label-based evals use `Mouse_Aging_subsetted_celltype_labelled.h5ad`.
- **Expected Output**: Ranked, uppercase, deduplicated gene list.
- **Ground Truth & Tolerance**: Lineage marker panels with Jaccard thresholds (≥ 0.50 for astrocyte/microglia, ≥ 0.38 for oligodendrocyte, ≥ 0.29 for ependymal/pericyte, ≥ 0.17 for endothelial) evaluated against the top-k predictions.
- **Core Test**: Confirms the agent can reproduce ageing signatures whether starting from scratch or from pre-labelled data.

**Lineage summaries**

- Astrocyte — [hierarchical](differential-expression/mouse_brain_merfish_astrocyte_age_up_onevsone_ttest_jaccard_v13_hier.json) / [label-based](differential-expression/mouse_brain_merfish_astrocyte_age_up_ttest_jaccard_v14_labels_de_only.json); reactive panel `{GFAP, VIM, C4B, C3, SERPINA3N, CXCL10, IL18, HIF3A}`; top 13; Jaccard ≥ 0.50.
- Endothelial — [hierarchical](differential-expression/mouse_brain_merfish_endothelial_age_up_onevsone_ttest_jaccard_v13_hier.json) / [label-based](differential-expression/mouse_brain_merfish_endothelial_age_up_ttest_jaccard_v14_labels_de_only.json); interferon sentinel `XDH`; top 6; Jaccard ≥ 0.17.
- Ependymal — [hierarchical](differential-expression/mouse_brain_merfish_ependymal_age_up_onevsone_ttest_jaccard_v13_hier.json) / [label-based](differential-expression/mouse_brain_merfish_ependymal_age_up_ttest_jaccard_v14_labels_de_only.json); interferon pair `{IFIT3, IFITM3}`; top 7; Jaccard ≥ 0.29.
- Microglia — [hierarchical](differential-expression/mouse_brain_merfish_microglia_age_up_onevsone_ttest_jaccard_v13_hier.json) / [label-based](differential-expression/mouse_brain_merfish_microglia_age_up_ttest_jaccard_v14_labels_de_only.json); DAA genes `{B2M, TREM2, CCL2, APOE, AXL, ITGAX, CD9, C1QA, C1QC, LYZ2, CTSS}`; top 16; Jaccard ≥ 0.50.
- Oligodendrocyte — [hierarchical](differential-expression/mouse_brain_merfish_oligodendrocyte_age_up_onevsone_ttest_jaccard_v13_hier.json) / [label-based](differential-expression/mouse_brain_merfish_oligodendrocyte_age_up_ttest_jaccard_v14_labels_de_only.json); inflammatory trio `{C4B, IL33, IL18}`; top 8; Jaccard ≥ 0.38.
- Pericyte — [hierarchical](differential-expression/mouse_brain_merfish_pericyte_age_up_onevsone_ttest_jaccard_v13_hier.json) / [label-based](differential-expression/mouse_brain_merfish_pericyte_age_up_ttest_jaccard_v14_labels_de_only.json); interferon pair `{IFIT3, IFITM3}`; top 7; Jaccard ≥ 0.29.

---

## Regional Composition Benchmarks (`vizgen/tissue_composition`)

- **What & Why**: Profiles must reproduce known cell-type mixtures for distinct anatomical regions, demonstrating that expression-only annotation yields anatomically faithful compositions without spatial priors.
- **Goal**: Using the shared hierarchy (Neuron vs Non-neuronal → lineage splits), compute 13 lineage percentages as `100 * count / total` for the specified tissue and emit `{"tissue": "...", "cell_type_percentages": {...}}`.
- **Input**: `latch://38438.account/Mouse_Aging_subsetted.h5ad`.
- **Expected Output**: JSON payload with the tissue name and raw (non-renormalised) percentages for the 13 required labels.
- **Ground Truth & Tolerance**: Region-specific reference percentages with absolute tolerance ±3.0 percentage points per lineage.
- **Core Test**: Ensures annotation, filtering, and counting preserve realistic regional compositions spanning neurons, glia, vasculature, and immune cells.

- [mouse_brain_pct_brain_ventricle_v4](tissue_composition/mouse_brain_pct_brain_ventricle_v4.json) — ventricular region; astrocyte ~24.65 %, MSN ~7.15 %.
- [mouse_brain_pct_corpus_callosum_v4](tissue_composition/mouse_brain_pct_corpus_callosum_v4.json) — corpus callosum; oligodendrocyte ~66.69 %.
- [mouse_brain_pct_layer_II_III_v4](tissue_composition/mouse_brain_pct_layer_II_III_v4.json) — cortical II/III; neuron ~55.45 %, inhibitory interneuron ~10.22 %.
- [mouse_brain_pct_layer_V_v4](tissue_composition/mouse_brain_pct_layer_V_v4.json) — cortical V; neuron ~48.54 %, inhibitory interneuron ~10.79 %.
- [mouse_brain_pct_layer_VI_v4](tissue_composition/mouse_brain_pct_layer_VI_v4.json) — cortical VI; neuron ~51.04 %, oligodendrocyte ~18.22 %.
- [mouse_brain_pct_olfactory_region_v4](tissue_composition/mouse_brain_pct_olfactory_region_v4.json) — olfactory region; neuron ~31.04 %, inhibitory interneuron ~23.85 %.
- [mouse_brain_pct_pia_mater_v4](tissue_composition/mouse_brain_pct_pia_mater_v4.json) — pia mater; astrocyte ~34.03 %, pericyte ~10.04 %, VLMC ~10.95 %.
- [mouse_brain_pct_striatum_v4](tissue_composition/mouse_brain_pct_striatum_v4.json) — striatum; medium spiny neuron ~64.60 %, oligodendrocyte ~10.70 %.

---

These evaluations collectively guarantee that MERFISH analyses of ageing mouse brain recover canonical cell-type markers, ageing-associated transcriptional programs, and anatomically realistic tissue compositions.

