# Cell Type Annotation Evaluation 

Accurate annotation requires verifying your results, not just assigning labels. This guide outlines the key metrics and visual checks needed to assess annotation quality and decide when to revise your approach.

- **IMPORTANT**: **Compute ALL metrics and visualization to self-evaluate the quality of cell type annotation.**

### 1. Cell Type Balance & Biological Plausibility

**Check:** 
- Do predicted proportions align with biological expectations of the dataset's **tissue**?
- Do **changes in proportions between conditions** look plausible biologically? 

### 2. Cluster-Celltype Alignment (Purity)

**Metric:** Mean cluster purity = average % of dominant cell type per cluster

**Calculation:** For each cluster, compute % belonging to most common cell type, then average across clusters.

**Benchmarks for ATAC-seq:**
- **>60%:** Excellent — predictions align well with clustering
- **50-60%:** Good — acceptable alignment for ATAC-seq
- **40-50%:** Moderate — some mismatch, but may be biological
- **<40%:** Poor — consider cluster-based annotation instead

**Why moderate purity is acceptable:** Gene activity (ATAC) captures different biology than what clustering may be based on (e.g., transcriptomics, spatial patterns). Moderate purity often reflects genuine biological heterogeneity within clusters.

### 3. Spatial Coherence

**Concept:** Similar cell types should be spatially co-localized in tissue sections. 

**Checks:**
- **Local homogeneity:** Do neighboring spots tend to have the same cell type?
- **Anatomical zones:** Do cell types localize to expected tissue regions? (e.g., hepatocytes in liver lobules, immune cells in portal triads) (if anatomical annotation is available in the AnnData object)

**Methods:**
- Calculate local neighborhood agreement (% of neighbors with same type)
- Check % anatomical zone overlap (if anatomical annotation is available)

**Red flags:**
- Unexpected mixing of incompatible cell types
- Predicted types in anatomically impossible locations


### 4. Validation of Discriminatory Marker per Cell Type

**Concept:** A valid discriminatory marker must be enriched in its target cell type compared to **EVERY other cell type individually**, not just show high average fold change.

**Why critical:** A marker can have high mean fold change while being **lower** in the target compared to one or more other cell types. Such markers fail to distinguish those specific pairs and are unsuitable for validation.

**Method:**

1. Calculate median expression per cell type for each marker
2. Compute fold change: target vs. each other cell type individually  
3. Identify **minimum fold change** (worst-case scenario)
4. Valid marker: minimum FC ≥1.0 (ideally ≥1.05)
5. For each predicted cell type, select its most discriminatory marker. Create a violin plot comparing its expression distribution across all cell types. Use one `w_plot` per cell type.

**Benchmarks for ATAC-seq:**

| Min FC Range | Assessment | Action |
|--------------|------------|--------|
| **>1.20×** | Excellent | Use confidently |
| **1.10-1.20×** | Good | Acceptable for ATAC-seq |
| **1.05-1.10×** | Moderate | Usable with caveats |
| **1.00-1.05×** | Poor | Note low confidence |
| **<1.00×** | Failed | Do not use; higher in other cell type(s) |

**Expected outcomes:**

- 70-80% of cell types: min FC >1.10×
- 10-20% of cell types: min FC 1.05-1.10× (marginal)
- 0-10% of cell types: no valid marker (annotation challenge)

**Red flags:**

- Immune cells often have lower min FC (1.03-1.15×) due to shared chromatin patterns
- Cell type with no valid marker indicates low annotation confidence
- Should flag for alternative methods (cluster-based DE, RNA integration)

**Rationale:**

Traditional validation using mean/median FC can mask critical failures where markers don't distinguish specific cell type pairs. Minimum FC ensures markers work in the **worst case**, essential for publication-quality validation and honest quality assessment.

### 5. Prediction Confidence Distribution

**Metric:** Mean confidence score per cell type (difference between top and second-best normalized score)

**Benchmarks:**
- **>0.4:** High confidence — clear cell type identity
- **0.25-0.4:** Good confidence — reliable predictions
- **0.15-0.25:** Moderate — acceptable but verify
- **<0.15:** Low — ambiguous predictions, use caution

**Cell type patterns:**
- Abundant types typically have higher confidence
- Rare types (<5% of cells) often have lower confidence
- Very low confidence may indicate marker insufficiency or mixed populations

**Action:** Flag cell types with mean confidence <0.15 as "low confidence" in downstream analyses.

### 6. Sample-Level Consistency

**For multi-sample datasets:** Check if cell type proportions are consistent or vary by condition.

**Expected patterns:**
- Control samples should have similar proportions
- Disease/treatment samples may differ (this is biology, not error)
- Technical replicates should be highly similar

**Red flag:** Wild variation between technical replicates suggests technical artifacts or insufficient normalization.

### 7. Reasonable Cell Type Proportion Changes Between Conditions

**For dataset with multiple conditions**: Check if cell type proportion changes _between_ conditions are biologically plausible. 

---

## Validation Summary Checklist

Before accepting predictions, verify:

- [ ] Cell type proportions are biologically plausible
- [ ] Mean cluster-celltype purity >50% (for ATAC-seq)
- [ ] Spatial patterns show coherent localization
- [ ] Markers show >1.2× enrichment in predicted types
- [ ] Mean confidence >0.20 overall
- [ ] No single cell type >80% unless expected
- [ ] Rare expected types are detected
- [ ] Sample-level consistency (for multi-sample data)
- [ ] Biologically reasonable cell type proportions across conditions

**If ≥5 of 8 criteria pass: Accept predictions**  
**If <5 criteria pass: Consider cluster-based annotation or threshold adjustment**
