This document provides ultimate guideline to load vizgen MERFISH datasets. Be honest and follow this document to the word.

## Experiment Setup**
- The user can either input a H5AD file or a directory containing the raw spatial images, 
- Check if users provide both an **H5AD file**   
    - If the user provided **only an H5AD file**, **ask** them for the spatial directory.  
    - The spatial directory should contains **pmtiles**, if so, dipslay the h5ad file and spatial_dir using the `w_h5` widget

- If the user provided raw images from the MERSCOPE instrument, 
    - Check that the images directory contain a `cell_metadata.csv`, `cell_by_gene.csv`, and `micron_to_mosaic_pixel_transform.csv`. Create an adata using, 
    ```python
    adata = sq.read.vizgen(
        path=dataset_dir,
        counts_file=counts_file,
        meta_file=meta_file,
        transformation_file=transformation_file,
    )
    ```
    - Use the `w_h5` widget to view the adata generated.  

- Spatial coordinates are typically found in `spatial` or `X_spatial`.

---

## Data Loading**
- Load the H5AD data using **Scanpy**.  
- Ensure the input is a valid **AnnData (H5AD)** object.  
- Always use **`w_h5`** that takes both **AnnData** and **spatial directory** as inputs.

---
