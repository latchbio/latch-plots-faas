<goal>
Perform differential abundance analysis to compare protein abundance between conditions (e.g., treated vs resting) within specific cell types or cell populations.
</goal>

<method>
1/ **Prepare Data Layer**
   - Ensure shifted_dsb layer exists (create if needed):
     ```python
     from pixelator.common.statistics import dsb_normalize
     import numpy as np
     import pandas as pd
     import scanpy as sc
     import seaborn as sns
     import matplotlib.pyplot as plt
     ```
     ```python
     adata.layers["dsb"] = dsb_normalize(
         adata.to_df(), 
         isotype_controls = isotype_controls
     )
     mins = adata.to_df("dsb").min()
     mins[mins > 0] = 0
     adata.layers["shifted_dsb"] = adata.to_df("dsb").sub(mins)
     ```

2/ **Subset Data and Define Conditions**
   - Allow user to select cell type or cell population to analyze (e.g., "CD8 T", "CD4 T", "B", etc.)
   - Subset data to selected cell type:
     ```python
     subset_adata = adata[adata.obs['cell_type'] == selected_cell_type, :]
     ```
   - Allow user to define condition groups based on sample names or metadata:
     - Specify condition names (e.g., "stimulated" vs "resting")
     - Define which samples belong to each condition
   - Assign condition labels:
     ```python
     subset_adata.obs['condition'] = 'condition_1'
     subset_adata.obs.loc[subset_adata.obs['sample'].str.contains('pattern'), 'condition'] = 'condition_2'
     ```

3/ **Differential Abundance Analysis**
   - Run rank genes groups analysis:
     ```python
     sc.tl.rank_genes_groups(
         subset_adata, 
         groupby='condition', 
         method='wilcoxon', 
         groups=['condition_2'], 
         reference='condition_1',
         layer='shifted_dsb'
     )
     ```
   - Extract differential expression results:
     ```python
     diff_exp_df = sc.get.rank_genes_groups_df(subset_adata, group=None)
     diff_exp_df["-log10(adjusted p-value)"] = -np.log10(diff_exp_df["pvals_adj"])
     ```
   - Allow user to specify p-value threshold (default: 0.01)
   - Mark significant genes:
     ```python
     p_value_threshold = 0.01
     diff_exp_df["Significant"] = diff_exp_df["pvals_adj"] < p_value_threshold
     ```
   - Display the differential abundance table

4/ **Volcano Plot Visualization**
   - Allow user to specify log fold change threshold for labeling (default: 1.5)
   - Generate volcano plot showing:
     - x-axis: log fold changes
     - y-axis: -log10(adjusted p-value)
     - Color coding based on significance
     - Labels for markers with absolute log fold change above threshold

5/ **Heatmap Visualization**
   - Allow user to specify number of top and bottom markers to display (default: 7 each)
   - Select top and bottom markers based on log fold changes
   - Calculate mean DSB values per sample for selected markers across all samples in the full dataset
   - Generate heatmap showing average DSB values per sample for the selected markers
</method>

<workflows>
</workflows>

<library>
</library>

<self_eval_criteria>
</self_eval_criteria>
