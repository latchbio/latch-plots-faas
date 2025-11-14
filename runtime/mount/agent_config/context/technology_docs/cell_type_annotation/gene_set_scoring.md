# Gene Set Scoring — Complete Workflow Guide

Gene set scoring is a **fast exploratory approach** for cell type annotation that evaluates how strongly each cell expresses curated marker gene panels. It works without requiring high-quality clustering and is ideal for initial annotation of spatial transcriptomics or ATAC-seq data.

**Key advantages:** Fast (5-10 min), no clustering dependency, interpretable  
**Key requirements:** Score normalization, balanced markers, proper thresholds

---

## Workflow Summary

1. **Select expected cell types** (5-10 major types for tissue)
2. **Curate 40-50 markers** per cell type from CellGuide
3. **Filter discriminatory markers** (median fold change threshold)
4. **Balance marker counts** across cell types
5. **⚠️ CRITICAL: Normalize scores** (z-score or min-max)
6. **Assign cell types** by highest normalized score
7. **⚠️ CRITICAL: Run the full evaluation suite to validate your own cell type assignments**: You must compute **every evaluation metric and visualization** defined in the evaluation section (proportions, cluster purity, spatial coherence, marker enrichment, confidence scores, sample consistency, and condition-wise proportion changes). **Do not skip or short-circuit any evaluation step.**
8. **Self-evaluate** based on metrics, reflect, and modify your work as needed. Be honest and objective.

---

## Step 1: Select Expected Cell Types

Choose **5-10 major cell types** realistic for your tissue and organism. Focus on broad categories rather than fine subtypes.

---

## Step 2: Curate Marker Panels from CellGuide

Extract **40-50 top-ranked markers** per cell type from the CellGuide database.

**Note**
- CellGuide names may not match the cell type names you’re looking for exactly.
- The database is comprehensive and includes broad categories, subtypes, and synonyms.
- Prioritize biological equivalence over exact string matches when searching for cell types and markers. 

**Why 40-50?**
- Expect 20-40% retention after discriminatory filtering
- Ensures ≥5 markers remain per cell type after filtering
- CellGuide ranks by marker strength — top markers are highest quality

**Data source:** `latch:///cellguide_marker_gene_database_per_celltype.json` from Latch Data

Filter to genes present in your dataset (`adata.var_names`)

---

## Step 3: Filter for Discriminatory Markers

Evaluate each marker's ability to distinguish cell types by computing **median fold change** across all cluster pairs.

### Recommended Thresholds by Assay Type

| Assay Type | Threshold | Expected Retention | Rationale |
|------------|-----------|-------------------|-----------|
| **ATAC-seq gene activity** | **1.5×** | 20-40% | Lower dynamic range, sparser signal |
| RNA-seq (scRNA/spatial) | 1.5-2.0× | 30-50% | Higher sensitivity and dynamic range |
| Protein (CITE-seq) | 2.0-3.0× | 40-60% | High specificity required |

**Key insight:** ATAC-seq gene activity scores are inherently sparser than RNA expression. The commonly used **1.5× threshold** is recommended but can be too stringent for some datasets. If that is the case, decrease the reshold **1.2× for ATAC-seq data.**

**Method:** For each marker, compute median expression per cluster, then calculate pairwise fold changes. Keep markers where median fold change across all pairs exceeds threshold.

**Quality check:** After filtering, verify **each cell type retains ≥3 markers AND there are ≥3 cell types**. If not, lower threshold to 1.2× or increase initial marker count to 50-60.

---

## Step 4: Balance Marker Counts (Target: 5 markers per cell type)

### ⚠️ CRITICAL: Prevent Scoring Bias

**Problem:** Unequal marker counts create systematic bias. A cell type with 20 markers will outscore one with 5 markers regardless of biological signal.

**Solution:** Equalize marker counts across all cell types.

**Method:**
1. Calculate target count (median or minimum of filtered marker counts)
2. Enforce minimum of 3 markers per cell type
3. Take top 5 markers by fold change for each cell type
4. Only include cell types with ≥3 final markers

**Result:** All cell types have equal "voting power" during scoring.

---

## Step 5: Compute Raw Gene Set Scores

Calculate **mean expression** across each cell type's marker panel for every cell.

**Formula:** For each cell and cell type, score = mean(expression of all markers for that cell type)

This produces a score matrix: cells × cell types

---

## Step 6: ⚠️ CRITICAL — Normalize Scores

### Why Normalization is Non-Negotiable

**Raw mean scores are NOT comparable across cell types** due to:
- Different baseline expression levels across marker sets
- Varying dynamic ranges of marker expression
- Technical biases in specific gene accessibility

**Without normalization:** Cell types with high baseline marker expression will dominate predictions (typically 80-99% mis-classification to a single type).

### Normalization Methods

**Z-score normalization (RECOMMENDED for ATAC-seq):**
- Centers each cell type at mean=0, std=1
- Formula: `(score - mean) / std` for each cell type
- Makes scores directly comparable across all cell types

**Min-max scaling (alternative):**
- Scales each cell type to 0-1 range
- Formula: `(score - min) / (max - min)`
- Useful when z-scores have outliers

**Rank-based (conservative):**
- Converts scores to percentile ranks
- Most robust but loses magnitude information

**For ATAC-seq data: Always use z-score normalization.**

---

## Step 7: Assign Cell Types

After normalization:
1. For each cell, identify the cell type with the **highest normalized score**
2. Calculate **prediction confidence** = difference between top and second-highest score
3. Add predictions and confidence to `adata.obs`

Higher confidence (>0.3) indicates clear cell type identity. Lower confidence (<0.15) suggests ambiguous or mixed populations.

---

## Step 8: Comprehensive Evaluation

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

### 4. Marker Gene Specificity

**Concept:** Markers should be enriched in their predicted cell types vs all others.

**Validation:** For each cell type and its top markers:
1. Calculate mean expression in predicted cells of that type
2. Calculate mean expression in all other cells
3. Compute fold change
4. Produce violin plots for users to visually inspect

**Benchmarks:**
- **>1.5×:** Excellent marker specificity
- **1.2-1.5×:** Good (expected for ATAC-seq)
- **<1.2×:** Poor — marker not truly discriminatory

**Check all cell types:** At least 60% of markers should show >1.2× enrichment

**Red flag:** If markers show similar or higher expression in other cell types, predictions are unreliable.

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

**If ≥4 of 8 criteria pass: Accept predictions**  
**If <4 criteria pass: Consider cluster-based annotation or threshold adjustment**

---

## Troubleshooting Common Issues

| Symptom | Likely Cause | Solution |
|---------|--------------|----------|
| 99% cells = one type | Missing normalization | Add z-score normalization (Step 6) |
| Only 2-3 types detected | Threshold too strict | Lower to 1.2× (ATAC) or 1.1× |
| Cluster purity <40% | Gene activity ≠ clustering | Use cluster-based approach |
| Rare types not detected | Insufficient markers | Increase initial markers to 50-60 |
| One type dominates (70-80%) | Marker imbalance | Balance marker counts (Step 4) |
| No spatial coherence | Poor quality predictions | Verify normalization and thresholds |
| Low marker specificity | Markers not discriminatory | Lower threshold or use different markers |
| Low overall confidence | Weak marker signal | Increase markers or try cluster-based |

---
