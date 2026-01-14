<goal>
Perform abundance normalization and denoising of PNA data using DSB (denoised and standardized by background) and CLR (Centered-Log-Ratio) methods, and visualize normalization results to assess transformation effects.
</goal>

<method>
1/ **Normalization Setup**
   - Import required Pixelator libraries:
     ```python
     from pixelator.mpx import read, simple_aggregate
     from pixelator.mpx.plot import density_scatter_plot
     from pixelator.common.statistics import dsb_normalize
     from pixelator.common.statistics import clr_transformation
     import numpy as np
     import pandas as pd
     import seaborn as sns
     import matplotlib.pyplot as plt
     ```
   - Get column names from adata: `columns = adata.to_df().columns`
   - Allow user to select isotype controls (default: ["mIgG1", "mIgG2a", "mIgG2b"])

2/ **Apply Normalization Methods**
   - Apply DSB normalization:
     ```python
     adata.layers["dsb"] = dsb_normalize(
         adata.to_df(), 
         isotype_controls = isotype_controls
     )
     ```
   - Create shifted DSB layer (shift minimum values to zero):
     ```python
     mins = adata.to_df("dsb").min()
     mins[mins > 0] = 0
     adata.layers["shifted_dsb"] = adata.to_df("dsb").sub(mins)
     ```
   - Apply CLR transformation:
     ```python
     adata.layers["clr"] = clr_transformation(
         adata.to_df(), axis=1
     )
     ```
   - Apply log1p transformation:
     ```python
     adata.layers['log1p'] = np.log1p(adata.to_df())
     ```

3/ **Normalization Visualization - Density Plots**
   - Allow user to select a marker of interest for evaluation (default: "CD19")
   - Store selected marker in variable `marker`
   - For log1p transformation:
     ```python
     layer_data_log = adata.to_df("log1p")
     data_log = layer_data_log.loc[:, marker]
     samples = adata.obs.loc[:, 'sample']
     plot_data_log = pd.DataFrame({marker : data_log, 'Samples' : samples})
     
     sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':2})
     pal = sns.color_palette("Set2", 12)
     g1 = sns.FacetGrid(plot_data_log, row="Samples", hue="Samples", aspect=9, height=1.2, palette=pal)
     g1.map(sns.kdeplot, marker, bw_adjust=.5, clip_on=False, fill=True, alpha=1, linewidth=1.5)
     
     def label(x, color, label):
         ax = plt.gca()
         ax.text(0, .2, label, fontweight="bold", color=color,
                 ha="left", va="center", transform=ax.transAxes)
     
     g1.map(label, marker)
     g1.figure.subplots_adjust(hspace=-.5)
     g1.set_titles("")
     g1.set(yticks=[], ylabel="")
     g1.despine(bottom=True, left=True)
     g1.fig.suptitle('Log1p Transformation')
     ```
     - Display the plot
   - For DSB transformation:
     ```python
     layer_data_dsb = adata.to_df("dsb")
     data_dsb = layer_data_dsb.loc[:, marker]
     plot_data_dsb = pd.DataFrame({marker : data_dsb, 'Samples' : samples})
     
     sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':2})
     pal = sns.color_palette("Set2", 12)
     g = sns.FacetGrid(plot_data_dsb, row="Samples", hue="Samples", aspect=9, height=1.2, palette=pal)
     g.map(sns.kdeplot, marker, bw_adjust=.5, clip_on=False, fill=True, alpha=1, linewidth=1.5)
     g.map(label, marker)
     g.figure.subplots_adjust(hspace=-.5)
     g.set_titles("")
     g.set(yticks=[], ylabel="")
     g.despine(bottom=True, left=True)
     g.fig.suptitle('dsb Transformation')
     ```
     - Display the plot

4/ **Density Scatter Plots**
   - Allow user to select two markers for density scatter plot comparison:
     - First marker (default: "CD19")
     - Second marker (default: "CD18")
   - Store selected markers in variables `marker_1` and `marker_2`
   - Generate density scatter plot for log1p layer:
     ```python
     fig_density_log_1p, ax_density_log_1p = density_scatter_plot(
         adata, 
         marker_1, 
         marker_2, 
         layer="log1p"
     )
     ```
   - Generate density scatter plot for DSB layer:
     ```python
     fig_density_dsb, ax_density_dsb = density_scatter_plot(
         adata, 
         marker_1, 
         marker_2, 
         layer="dsb"
     )
     ```
   - Display both plots side by side with labels "log1p" and "dsb"
</method>

<workflows>
</workflows>

<library>
</library>

<self_eval_criteria>
</self_eval_criteria>
