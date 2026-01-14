<goal>
Perform quality control and data integration analysis on Pixelator data, including cell calling validation, sample-level cell counts, QC metric distributions, and antibody count distribution outlier detection.
</goal>

<method>
1/ **Edgerank Plot - Cell Calling Quality Control**
   - Import Pixelator plotting function:
     ```python
     from pixelator.pna.plot import molecule_rank_plot
     import seaborn as sns
     ```
   - Set seaborn style: `sns.set_style("whitegrid")`
   - Create an interactive slider for molecule filter threshold:
     - Default value: 10000
     - Min: 0
     - Max: 50000
     - Step: 100
   - Create `molecule_rank_df` by copying `adata.obs[["condition", "n_umi"]]`
   - Compute rank column: `molecule_rank_df.groupby(["condition"])["n_umi"].rank(ascending=False, method="first")`
   - Generate the molecule rank plot using `molecule_rank_plot(molecule_rank_df, group_by="condition")`
   - Add a horizontal reference line at the threshold value: `ax_edge_rank.axhline(int(significance_threshold.value), linestyle="--")`
   - Display the plot with label "Cell calling: Edge rank plot"

2/ **Cell Counts per Sample**
   - Import required libraries:
     ```python
     import plotly.express as px
     ```
   - Create `cells_per_sample_df` by grouping `adata.obs` by "sample" and computing size:
     ```python
     cells_per_sample_df = (
         adata.obs.groupby("sample").size().to_frame(name="size").reset_index()
     )
     ```
   - Print unique samples: `print(adata.obs['sample'].unique())`
   - Create a bar plot using `px.bar`:
     - Data: `cells_per_sample_df`
     - x-axis: "sample"
     - y-axis: "size"
     - color: "sample"
     - title: "Counts per Sample"
   - Format the plot:
     - Update traces with `textposition="outside"` and `width=0.4`
     - Update layout with:
       - xaxis_title: "Sample"
       - yaxis_title: "Count"
       - showlegend: False
   - Display the plot with label "No. of cells per sample"
   - Display the data table with `cells_per_sample_df`

3/ **QC Distribution Plots**
   - Import required libraries:
     ```python
     import seaborn as sns
     ```
   - Compute `reads_per_umi` metric: `adata.obs["reads_per_umi"] = adata.obs["reads_in_component"]/adata.obs["n_umi"]`
   - Create `metrics_per_sample_df` by melting the observation dataframe:
     ```python
     metrics_per_sample_df = adata.obs[
         ["condition", "n_umi", "reads_per_umi", "average_k_core"]
     ].melt(id_vars=["condition"])
     ```
   - Generate violin plots using `sns.catplot`:
     - data: `metrics_per_sample_df`
     - x: "condition"
     - col: "variable"
     - y: "value"
     - kind: "violin"
     - sharex: True
     - sharey: False
     - margin_titles: True
     - hue: "condition"
   - Set column titles template: `.set_titles(col_template="{col_name}")`
   - Display the plot with label "QC Plots"

4/ **Antibody Count Distribution Outlier Removal**
   - Import required libraries:
     ```python
     import numpy as np
     import seaborn as sns
     ```
   - Compute log10-transformed n_umi: `adata.obs["n_umi_log10"] = np.log10(adata.obs["n_umi"])`
   - Create `tau_metrics_df` by selecting relevant columns:
     ```python
     tau_metrics_df = adata.obs[["condition", "tau", "n_umi_log10", "tau_type"]]
     ```
   - Generate scatter plot using `sns.relplot`:
     - data: `tau_metrics_df`
     - x: "tau"
     - y: "n_umi_log10"
     - col: "condition"
     - kind: "scatter"
     - hue: "tau_type"
   - Display the plot with label "Antibody count distribution outlier removal"
</method>

<workflows>
</workflows>

<library>
</library>

<self_eval_criteria>
</self_eval_criteria>
