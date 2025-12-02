## Feature Selection

— identify **highly variable genes (HVGs)**.
- Use sc.pp.highly_variable_genes(..., subset=True) to subset the data to the selected HVGs.
- Allow the user to specify the number of HVGs to retain (default: 2,000 genes).
- Always generate diagnostic plots using Scanpy for visualization.
