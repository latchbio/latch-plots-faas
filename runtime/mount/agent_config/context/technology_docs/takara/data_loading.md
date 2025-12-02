## Data Loading

- Use your knowledge of the kit type - Seeker 3x3, Seeker 10x10, or Trekker - to modify the code used for background removal. Ask for this information if you don't already have it.
- Ask whether their H5AD contains **one** or **multiple samples**. Spatial coordinates are typically found in **spatial** or **X_spatial**
- Load data using **Scanpy**.  Ensure the input is a valid **AnnData (H5AD)** object.  **ALWAYS** use ```w_h5``` for displaying anndata
