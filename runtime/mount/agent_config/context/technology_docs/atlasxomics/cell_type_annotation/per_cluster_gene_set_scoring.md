## Cluster-based Gene Set Scoring

### Overview
This protocol identifies cell types in AtlasXOmics spatial ATAC-seq data using **literature-based marker genes**. AtlasXOmics measures chromatin accessibility and imputes gene activity scores, not direct RNA expression.

--

## Workflow Summary

1. **Define expected cell types** — Collect literature-based marker genes (≥10 per type) suitable for ATAC-derived gene activity.

2. **Check marker availability** — Verify markers exist in `adata.var_names`; update lists if coverage is low.

3. **Score clusters** — Compute mean gene activity and Z-scores for each marker set and compare Z-scores **across clusters** to identify enriched cell types.

4. **Assign annotations** — Label each cluster with the cell type showing the highest relative marker score.

5. **Add results to AnnData** — Write primary cell type assignments (and optional detailed labels) to `adata.obs`.

6. **⚠️ CRITICAL: Evaluate your work** — Read `technology_docs/atlasxomics/cell_type_annotation/evals.md`. Compute **all** required metrics: proportions, purity, spatial coherence, marker enrichment, confidence, sample consistency, condition effects.  

7. **Revise** — If results are weak, adjust markers, thresholds, or scoring logic, or switch to another annotation approach.

---

### Key AtlasXOmics Considerations

**Data Type:**
- **Gene activity scores** derived from ATAC-seq peaks mapped to gene bodies
- **NOT direct RNA expression** - accessibility proxy for gene activity
- Scores are typically **5-10% lower** than RNA-seq equivalents
- Lowly expressed genes may not show accessibility signal

**Expected h5ad Structure:**
- `adata.X`: Gene activity matrix (imputed from chromatin accessibility)
- `adata.obs['cluster']`: Cluster assignments from upstream analysis
- `adata.obs['n_fragment']`: Total ATAC fragments per cell/spot
- `adata.obs['tsse']`: TSS enrichment score (quality metric)
- `adata.obsm['X_umap']`: UMAP coordinates
- `adata.obsm['spatial']`: Spatial coordinates (row, col or xcor, ycor)

---

### Annotation Workflow

#### **1. Identify Expected Cell Types**

Research the major cell types for your tissue and their canonical markers. Use 10+ markers per cell type for robust scoring.

**Marker Sources:**
- CellGuide database: `cellguide_marker_gene_database_per_celltype.json` on Latch
- PubMed: Search "[tissue] cell type markers"
- Recent single-cell atlases for your tissue
- Key review papers

**Strategy:** Focus on genes with well-characterized regulatory elements that will show strong ATAC signal.

#### **2. Check Marker Presence**

Verify which of your markers exist in `adata.var_names`. If <5 markers present for a cell type, find additional markers or reconsider that cell type.

#### **3. Calculate Expression Scores (Z-score for Relative Enrichment)**

For each cluster, compute the mean marker expression for each cell type. This produces a score matrix where rows represent clusters and columns represent cell types.

**Key principle:** Cell type identity depends on **relative enrichment across clusters**, not absolute magnitude. Some marker sets are globally high everywhere, which can cause incorrect assignments if raw means are used.

**Required Approach**

- Compute mean marker expression per cluster.  
- **Z-score normalize each cell type across clusters** to convert raw scores into enrichment values.  
- Assign the cell type whose marker set shows the **highest Z-score** in that cluster.

**Why Z-scores Matter**

- Highlights where a marker set is **uniquely enriched**, rather than simply abundant.  
- Corrects for marker sets that are high across all clusters.  
- Produces biologically meaningful, relative signals that reflect true cell type specificity.

#### **4. Assign Primary Cell Types**

Assign each cluster to the cell type with the highest marker score. Clusters with similar low scores across all types need investigation.

#### **5. Investigate Unknown Clusters**

**Check quality first:**
- Low fragment count and TSS enrichment? → Technical failure, flag to user. 
- Good quality but low scores? → Proceed to investigation

**Find identifying features:**
- Extract top 20 highly expressed genes
- Check if they suggest a cell type you didn't test

**Test additional categories if needed:**
- Immune cells: Cd3d/e/g (T cells), Cd19/Ms4a1 (B cells), Cd68/Adgre1 (macrophages), Ncr1/Klrb1c (NK)
- Stromal: Col1a1/Col1a2 (fibroblasts), Rgs5/Pdgfrb (pericytes), Acta2/Tagln (smooth muscle)
- Endothelial: Pecam1, Cdh5, Kdr (present in most tissues)

**Quiescent/low activity states:**
- Good quality metrics but uniformly 5-10% lower expression across ALL markers
- Not a technical artifact - biological phenomenon
- Likely G0 cell cycle phase or metabolically inactive state

#### **6. Subtype Analysis (Optional)**

For major cell types spanning multiple clusters, define subtype markers and re-score.

**Example - Hepatocyte zonation:**
- Periportal (Zone 1): Hal, Ass1, Arg1, Pck1 (urea cycle, gluconeogenesis)
- Mid-zonal (Zone 2): Hsd17b13, Cyp2c29
- Pericentral (Zone 3): Glul, Cyp2e1, Cyp1a2 (glutamine synthesis, detox)

Calculate pericentral/periportal ratio for each cluster to assign zones.

#### **7. Spatial Validation**

Use `w_h5` interactive viewer to visualize:
- Cell type distribution on spatial coordinates
- Individual marker gene expression patterns
- Co-localization with anatomical features
- Spatial gradients (e.g., zonation patterns)

**Validation questions:**
- Do cell types cluster appropriately? (e.g., endothelial outlining vessels)
- Do markers show expected spatial patterns?
- Are there anatomical landmarks matching your assignments?

---

### AtlasXOmics-Specific Tips

**Why ATAC-seq changes things:**
- **Transcription factors** and **chromatin remodelers** often show strong signal
- **Structural proteins** with constitutive expression may be weaker
- **Secreted proteins** may be less detectable than in RNA-seq
- Focus on markers with **strong enhancer/promoter activity**

**Best marker types for ATAC:**
- Cell identity transcription factors (e.g., Hnf4a for hepatocytes, Gfap for astrocytes)
- Surface receptors with active transcription (e.g., Kdr, Pecam1)
- Metabolic enzymes with regulated expression (e.g., Cyp2e1, Glul)

**Avoid relying solely on:**
- Highly secreted proteins (may not correlate with accessibility)
- Constitutively expressed housekeeping genes
- Post-transcriptionally regulated genes

---

### Common Pitfalls

**1. Absolute vs Relative Scoring**
Don't expect high absolute scores like RNA-seq. Compare across clusters to find enrichment.

**2. Mixed Populations**
AtlasXOmics tiles are ~5-10μm, may contain parts of multiple cells. Mixed scores are possible in transition zones.

**3. Rare Cell Types**
May form small clusters with weak signals. Check spatial distribution for anatomical clues.

**4. Doublets**
High fragment count + mixed markers → potential doublet. Check if spatial distribution suggests technical artifact.

---

### Output Annotations

Add to `adata.obs`:
- `cell_type`: Primary cell type assignment
- `cell_type_detailed`: Optional subtype information (e.g., "Hepatocytes (Zone 1)")
- Document marker lists and rationale for reproducibility

---

### Success Criteria

✓ All major tissue cell types identified with clear marker support  
✓ Unknown clusters explained (low quality, quiescent state, or rare type)  
✓ Spatial patterns validated and biologically plausible  
✓ Annotations reproducible with documented marker sources  
✓ Quality metrics checked for all assignments

---

### Tissue-Specific Notes

**Liver:** Expect hepatocytes (50-70%), endothelial/LSEC (5-15%), Kupffer (5-10%), HSC (3-8%), cholangiocytes (1-3%). Hepatocytes show zonation - test for it.

**Brain:** Neurons dominant (40-60%), astrocytes (10-20%), oligodendrocytes (5-15%), microglia (5-10%). Watch for regional heterogeneity.

**Kidney:** Tubular epithelial cells dominant (60-80%), glomerular cells (5-10%), immune infiltrate variable. Segment-specific markers crucial.

**Lung:** Epithelial cells (40-60%, multiple types), endothelial (15-25%), immune (10-20%). AT1 vs AT2 distinction important.

**Adapt marker lists to your specific tissue biology.**
