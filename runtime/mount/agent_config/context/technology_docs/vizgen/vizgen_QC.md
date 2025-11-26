This document provides ultimate guideline to perform quality control analysis of vizgen MERFISH datasets. Be honest and follow this document to the word.

## Quality Control & Filtering**

- Use **Scanpy** to compute QC metrics.  
- **ALWAYS** include standard metrics such as:
  - Total counts per cell  
  - Number of detected genes per cell  
  - Mitochondrial gene percentage  
  - **Cell volume distribution** (derived from polygon area or segmentation mask; visualize histogram). You can never ignore this plot. 
- **ALWAYS** make **histograms** of QC metrics.  
- **ALWAYS** expose filtering parameters using `w_text_input`.  
- **ALWAYS** confirm with the user before applying filters using **lplots widgets**.  
- **ALWAYS** report how many cells will be removed before filtering.  
- If the user declines, **revert to the original AnnData**.
