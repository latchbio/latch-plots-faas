This document provides ultimate guideline to perform spatial analysis on vizgen MERFISH datasets. Be honest and follow this document to the word.

All spatial analyses must use **Squidpy** functions (`sq.gr`, `sq.pl`) and rely on spatial coordinates within the AnnData object.

### ** Centrality Analysis**
- Compute **spatial centrality metrics** (e.g., degree, betweenness, closeness) using `sq.gr.centrality_scores`. If these centrality scores do not exist, compute these sentrality scores and plot them.
- Plot **centrality plots** overlayed on spatial coordinates.  
- Allow the user to select which metric(s) to visualize. 
```python
sq.gr.centrality_scores(adata, cluster_key="leiden")
```

### ** Co-occurrence Analysis**
- Compute **pairwise co-occurrence scores** between clusters using `sq.gr.co_occurrence`.  
- **Allow users to specify** which clusters or cell types to compare using **lplots widgets**.  
- Display **heatmaps** or **bar plots** of co-occurrence frequencies using `sq.pl.co_occurrence`.
- Compute the cooccurrence similar to the snippet below and use it generate a new figure. 
```
sq.gr.co_occurrence(adata, cluster_key="leiden", n_jobs = 4)

# Extract co-occurrence results
# The result is 3D: (n_clusters, n_clusters, n_intervals)
# We'll average across intervals to get a single 2D matrix
cooc_data = adata.uns['leiden_co_occurrence']['occ']
cooc_mean = np.mean(cooc_data, axis=2)  # Average across spatial intervals

# Create heatmap manually using seaborn
fig_cooc, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(
    cooc_mean,
    cmap='RdBu_r',
    center=0,
    annot=False,
    square=True,
    cbar_kws={'label': 'Mean Co-occurrence Score'},
    ax=ax
)
```
### ** Neighborhood Enrichment**
- Compute **neighborhood enrichment** between clusters using `sq.gr.nhood_enrichment`.  
- Visualize using **heatmaps** and **spatial overlays** (`sq.pl.nhood_enrichment`).  
- Optionally display **z-score–normalized enrichment**.

### ** Ripley’s Statistics**
- Compute **Ripley’s L** statistics using `sq.gr.ripley`.  
- Plot **observed vs. expected K(r)** curves for major clusters.  
- Highlight clusters showing **significant spatial aggregation**.
- Plot the Ripley's statistic for **EVERY** cluster using plotly. 
- Allow the user to select between L,F anf G statistics. Default (F)
```python
mode = "F"
sq.gr.ripley(adata, cluster_key="leiden", mode=mode)

# Access the ripley results - key format is 'leiden_ripley_L'
ripley_key = f"leiden_ripley_{mode}"
if ripley_key not in adata.uns:
    # Try alternative key format
    ripley_key = "leiden_ripley"

ripley_results = adata.uns[ripley_key]
n_clusters = adata.obs['leiden'].nunique()

# Create line plot for all clusters
df = ripley_results[f"{mode}_stat"]   # tidy dataframe with: bins, leiden, stats

# Create Plotly figure
fig_ripley = go.Figure()

# Iterate over clusters present in the DF
for cluster_id in sorted(df["leiden"].unique()):
    sub = df[df["leiden"] == cluster_id]
    
    fig_ripley.add_trace(go.Scatter(
        x=sub["bins"],
        y=sub["stats"],
        mode='lines',
        name=f'Cluster {cluster_id}',
        line=dict(width=2)
    ))

# Add reference 0-line (CSR expectation for L-statistic)
fig_ripley.add_hline(
    y=0,
    line_dash="dash",
    line_color="black",
    annotation_text="Random Distribution"
)

fig_ripley.update_layout(
    title="Ripley's L Statistic for All Clusters",
    xaxis_title="Spatial Distance (r)",
    yaxis_title=f"Ripley's {mode}(r)",
    height=600,
    showlegend=True,
    legend=dict(
        orientation="v",
        yanchor="top",
        y=1,
        xanchor="left",
        x=1.02
    )
)
```

### ** Spatial Autocorrelation (Moran’s I)**
- Compute **Moran’s I** for gene expression spatial autocorrelation using `sq.gr.spatial_autocorr`.  
- Rank genes by Moran’s I score.  
- Display a **table of top 10 spatially autocorrelated genes**.  
- Optionally, show **spatial expression plots** for these genes with `sq.pl.spatial_scatter`.
