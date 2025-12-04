# Single-Cell Reference Curation

**Task:** Curate a single-cell AnnData reference from a paper.  
**Requirement:** Keep **raw counts in `adata.X`** for cell2location compatibility.  
**Global Rule:** At every step, the agent must output tables, plots (when applicable), and a brief markdown summary.

---

## Behavior Modes

- **Mode Check (ALWAYS FIRST):**  
  Check the current behavior mode is `"<turn_structure></turn_structure>"` and confirm this internally before continuing.

- **Step-by-step:**  
  After *each step*, STOP and ask the user for feedback before proceeding.

- **Proactive / Default:**  
  Execute the entire workflow end-to-end without pauses unless a critical error occurs.

---

## 1. Paper Input
Ask the user to upload a **text file** or paste the paper text. PDFs are not accepted.

## 2. Extract Cell Types
Parse the text to extract **author-defined cell types** and any marker/cluster descriptions.

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
Inspect top markers per cluster and assign cell types consistent with the paper.

## 7. Marker Gene Validation
Validate annotations using canonical marker genes (e.g., violin plots).  
Assess whether markers cleanly separate expected cell types.

## 8. Self-Evaluation
Reflect on whether clusters and markers match expected biology.  
If validation fails, clustering is insufficient.

## 9. Iterate
Re-run clustering with different resolutions (e.g., 0.5, 1.5, 2.0).  
Re-run DE → annotation → validation.  
Repeat until annotations are **biologically sound**.

## 10. Finalize Reference
Save AnnData with:  
- Raw counts preserved in `.X`  
- Final cell type annotations  
- Clean metadata.
