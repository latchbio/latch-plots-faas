# Single-Cell Reference Curation

**Task:** Curate a single-cell AnnData reference from a paper.  
**Requirement:** Keep **raw counts in `adata.X`** for cell2location compatibility.  
**Global Rule:** At every step, the agent must output tables, plots (when applicable), and a brief markdown summary.

---

## Behavior Modes

- **Mode Check (ALWAYS FIRST):**  
  Check the current behavior mode is `"<turn_structure></turn_structure>"` and confirm this internally before continuing.

- **Step-by-step:**  
  After *each step*, the agent must **STOP** and explicitly ask the user for feedback or approval before continuing.

- **Proactive / Default:**  
  First, **review the full plan** and request **all missing inputs or clarifications upfront**.  
  Once inputs are confirmed, execute the **entire workflow end-to-end** without pausing unless a critical error occurs.

---

## 1. Paper Input
Request that the user upload a **text file** or paste the paper text.  
(**PDFs are not accepted.**)

## 2. Extract Cell Types
Parse the text to extract **author-defined cell types** and any described markers or clusters.  
Ask the user to confirm whether these are the exact cell types they want annotated in the curated reference.

## 3. Preprocessing
Load raw data.  
Apply QC filters.  
Normalize **only for clustering** (retain raw `.X`).  
Compute HVGs → PCA → neighbor graph.

## 4. Clustering
Run Leiden clustering (initial resolution = **1.0**).

## 5. Differential Expression
Perform **1-vs-all DE** using Scanpy:  
`method="t-test_overestim_var"`.

## 6. First-Round Annotation

Use **Thinking mode** and **manually inspect the top 10–20 markers** for each cluster to assign cell types.

**Never** rely on heuristic code or automated “top marker” rules — these miss nuance and routinely mislabel spatial datasets. **Manual reasoning is mandatory**.

## 7. Marker Gene Validation (Quantitative)

Compute marker expression **in code** and judge quality using numbers:

For each marker × cell type:
- Mean/median expression per cluster  
- Fraction of cells > 0  
- Relative enrichment (fold-change or z-score) vs all other clusters  

Validation passes **only if** the expected cluster has **clear, strong enrichment** and off-target clusters do not.

## 8. Critical Self-Evaluation (Be Skeptical)
Your default stance is **critical, not confirmatory**.

For each annotated cell type:
- **Evidence FOR:** List specific markers + quantitative stats.  
- **Evidence AGAINST:** List any conflicting markers, weak signals, or multi-cluster expression.  
- **Verdict:**  
  - **Pass** only if evidence is strong and unambiguous.  
  - Otherwise mark **Needs Revision** and redo clustering/annotation.

Never accept your own result without trying to **disprove it first**.

## 9. Iterate
- Examine all cluster annotations for any potential mislabels. 
- Consider re-running clustering with different resolutions (e.g., 0.5, 1.5, 2.0). Re-run step 5 → step 8. Follow every rule at each step; do not skip.
- Repeat this process (3+ times) until annotations are **biologically sound**.

## 10. Finalize Reference
Save AnnData with:  
- Raw counts preserved in `.X`  
- Final cell type annotations  
