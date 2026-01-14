<goal>
Perform differential expression analysis to identify marker genes/proteins for each cluster and assign cell-type labels using both heatmap-based analysis and marker gating approaches, then compare the two methods.
</goal>

<method>
1/ **Prepare Data Layer**
   - Create shifted DSB layer if not already present:
     ```python
     mins = adata.to_df("dsb").min()
     mins[mins > 0] = 0
     adata.layers["shifted_dsb"] = adata.to_df("dsb").sub(mins)
     ```

2/ **Differential Expression Analysis**
   - Import required libraries:
     ```python
     import scanpy as sc
     import numpy as np
     import pandas as pd
     import seaborn as sns
     ```
   - Run rank genes groups analysis:
     ```python
     sc.tl.rank_genes_groups(
         adata, "leiden", method="wilcoxon", layer="shifted_dsb"
     )
     ```
   - Extract differential expression results:
     ```python
     diff_exp_df = sc.get.rank_genes_groups_df(adata, group=None)
     diff_exp_df["-log10(adjusted p-value)"] = -np.log10(diff_exp_df["pvals_adj"])
     ```

3/ **Significance Filtering**
   - Allow user to specify p-value threshold
   - Mark significant genes:
     ```python
     p_value_threshold = float(p_val_threshold)
     diff_exp_df["Significant"] = diff_exp_df["pvals_adj"] < p_value_threshold
     ```
   - Display the differential expression table

4/ **Volcano Plot Visualization**
   - Allow user to select a cluster ID
   - Generate volcano plot for the selected cluster showing:
     - x-axis: log fold changes
     - y-axis: -log10(adjusted p-value)
     - Color coding based on significance threshold

5/ **Heatmap Generation**
   - Create pivot table for heatmap:
     ```python
     df = diff_exp_df.pivot(index=["names"], columns=["group"], values=["logfoldchanges"])
     ```
   - Filter markers for heatmap (log2 fold change > 1 and significant in at least one cluster):
     ```python
     markers_for_heatmap = set(
         diff_exp_df[
             (np.abs(diff_exp_df["logfoldchanges"]) > 1) & diff_exp_df["Significant"]
         ]["names"]
     )
     df = df[df.index.isin(markers_for_heatmap)]
     ```
   - Format column names:
     ```python
     df.columns = [cluster for _, cluster in df.columns]
     ```
   - Generate clustermap:
     ```python
     fig = sns.clustermap(df, yticklabels=True, xticklabels=True, linewidths=0.1, cmap="vlag")
     ```
   - Display the heatmap

6/ **Cell-Type Assignment - Heatmap-Based**
   - Analyze the heatmap to identify enrichment patterns:
     - Clusters are displayed along the x-axis
     - Proteins/markers are displayed along the y-axis
     - Positive values indicate enrichment, negative values indicate depletion
   - Use the enrichment and depletion patterns of known cell-type markers to assign cell-type labels to clusters
   - Store assignments in `adata.obs["cell_type_heatmap"]` or similar annotation column

7/ **Marker-Based Gating for Cell-Type Assignment**
   - Define cell-type marker sets:
     ```python
     celltype_marker_sets = {
         "CD8T": ["CD3e", "CD8", "-CD19", "-CD4"],
         "CD4T": ["CD3e", "CD4", "-CD19", "-CD8"],
         "B": ["CD19", "CD20", "-CD3e", "-CD11c"],
         "Monocyte": ["CD11c", "CD163", "-CD3e", "-CD19"],
         "NK": ["CD335", "CD56", "-CD3e", "-CD19"],
         "Platelet": ["CD36", "CD41", "CD62P", "-HLA-DR"],
     }
     ```
   - Define gate calculation function:
     ```python
     def calc_gates_passed(expression_data: pd.DataFrame, marker_sets):
         gate_threshold = 1.5
         celltypes = list(marker_sets.keys())
         celltype_score = pd.DataFrame(index=expression_data.index, columns = celltypes)
         for ct in celltypes:
             celltype_score.loc[:, ct] = 0
             for marker in marker_sets[ct]:
                 if marker[0] == '-': # Marker should be missing from celltype ct
                     celltype_score.loc[:, ct] += (expression_data.loc[:, marker[1:]] < gate_threshold)
                 else:
                     celltype_score.loc[:, ct] += (expression_data.loc[:, marker] > gate_threshold)
             celltype_score.loc[:, ct] = celltype_score.loc[:, ct] / len(marker_sets[ct])
         return celltype_score
     ```
   - Allow user to specify required gate pass ratio (default: 0.75)
   - Calculate cell-type scores and assign labels:
     ```python
     required_gate_pass_ratio = 0.75
     celltype_scores = calc_gates_passed(adata.to_df("dsb"), celltype_marker_sets)
     adata.obs["cell_type_gate"] = celltype_scores.idxmax(axis=1).values
     adata.obs.loc[celltype_scores.max(axis=1) < required_gate_pass_ratio, "cell_type_gate"] = "other"
     ```

8/ **Comparison of Cell-Type Assignment Methods**
   - Compare annotations between heatmap-based and marker gating methods
   - Generate and display a confusion matrix showing agreement between the two assignment methods
   - Visualize the confusion matrix to assess concordance
</method>

<workflows>
</workflows>

<library>
</library>

<self_eval_criteria>
</self_eval_criteria>
