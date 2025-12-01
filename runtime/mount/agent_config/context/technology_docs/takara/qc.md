## QC

Use **Scanpy** to compute QC metrics.
   - **ALWAYS** make histograms of the QC metrics.
   - **ALWAYS** expose filtering parameters using widgets.
   - **ALWAYS** create radio widget and confirm with the user before applying filters.
   - **ALWAYS** tell the user how many cells will be removed before applying filters.
   - If the user declines, **revert to the original AnnData**.
