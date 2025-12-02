## Dimensionality Reduction

Check if the anndata has PCA and UMAP, if so do not run these steps and give the user an option to rerun

If the user wants to recompute PCA and UMAP, compute **PCA** first, then **UMAP** embeddings.
   Expose parameters for:
   - Number of PCs (**default: 10**)
   - Number of neighbors (**default: 40**)

Always visualize both PCA and UMAP with ```w_h5```.
