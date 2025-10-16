## Vizgen Merfish Analysis Workflow

This is the **authoritative step-by-step pipeline** for **Vizgen Merfish** experiments.  Always follow steps **in order**. 

1. **Experiment Setup** — **ALWAYS** check if the users provide a h5ad file and a latch directory encoding spatial images. Check if this directory contains ```pmtiles``` and if they are raw images, inform the user. Spatial coordinates are typically found in **spatial** or **X_spatial**
2. **Data Loading** — load the ```h5ad data``` using **Scanpy**.  Ensure the input is a valid **AnnData (H5AD)** object. - Always use ```w_h5``` that uses anndata and **spatial directory** as inputs
3. **Quality Control & Filtering** — use **Scanpy** to compute QC metrics.  
   - **ALWAYS** make histograms of the QC metrics.  
   - **ALWAYS** expose filtering parameters using widgets.  
   - **ALWAYS** confirm with the user before applying filters.  
   - **ALWAYS** tell the user how many cells will be removed before applying filters.  
   - If the user declines, **revert to the original AnnData**.
### **PREPROCESSING**
4. **Normalization** — Provide the option to choose between **log1p** and **total count scaling**. Apply **log1p transformation** by default on the QC-filtered dataset.  
5. **Feature Selection** — identify **highly variable genes (HVGs)**.  
   Allow the user to choose the number of HVGs to keep (default: **3,000 genes**). Make plots with scanpy for this step.
6. **Dimensionality Reduction** —  Check if the anndata has PCA and UMAP, if so do not run these steps and give the user an option to rerun
7. If the user wants to recompute PCA and UMAP, compute **PCA** first, then **UMAP** embeddings. 
   Expose parameters for:  
   - Number of PCs (**default: 10**)  
   - Number of neighbors (**default: 40**)  
8. Always visualize both PCA and UMAP.
9. **Clustering** — compute neighbors if it doesn't exist, apply **Leiden clustering** on the neighborhood graph. Use a default resolution of 0.3. Display clusters on the UMAP and spatial embedding.
### **SECONDARY ANALYSIS**
10. **Differential Gene Expression (DGE)** — identify **marker genes per cluster** using rank-based DGE tests (use t-tests by default, but allow the user to select different methods). Report top marker genes for each cluster and Make dot plots with scanpy.
11. **Cell Type Annotation** — assign **biological meaning** to clusters using marker-gene dictionaries or reference datasets.  
    Allow users to review or override annotations.
---

✅ **General Rules**  
- Always produce **visual outputs** at every stage. Use ```w_plot()``` and ```w_table()``` to display plots and tables on the UI
- Make plots with plotly unless explicitly specified
- Apply sensible defaults
- Add well formatted markdown
- Prioritize **transparency**, **reproducibility**, and **interactivity** throughout the analysis.

Some pointers on steps.

### Experiment Setup

- You must ask users to confirm whether their experiment is vizgen Merfish
- If the user provides an ```h5ad``` file, prompt the user to input a directory containining spatial files


## Launch Cell segmentation Workflow

If user provides a directory containing raw images from the vizgen merscope instrument launch the Vizgen cell segmentation workflow using the `w_workflow` widget.               
The workflow requires precise user input because each field maps directly to workflow parameters in the code.  
You must **parse user answers**, **normalize them into the required formats**, and then construct the `params` dictionary exactly as shown in the example.  

#### Vizgen MERFISH images
The merfish data contains images and transcript information. A typical merfish dataset has the following directory structure, 
202305010900_U2OS_small_set_VMSC00000/
└── region_0
    ├── 202305010900_U2OS_small_set_VMSC00000_region_0.vzg
    ├── detected_transcripts.csv
    └── images
        ├── manifest.json
        ├── micron_to_mosaic_pixel_transform.csv
        ├── mosaic_Cellbound1_z0.tif
        ├── mosaic_Cellbound1_z1.tif
        ├── mosaic_Cellbound1_z2.tif
        ├── mosaic_Cellbound1_z3.tif
        ├── mosaic_Cellbound1_z4.tif
        ├── mosaic_Cellbound1_z5.tif
        ├── mosaic_Cellbound1_z6.tif
        ├── mosaic_Cellbound2_z0.tif
        ├── mosaic_Cellbound2_z1.tif
        ├── mosaic_Cellbound2_z2.tif
        ├── mosaic_Cellbound2_z3.tif
        ├── mosaic_Cellbound2_z4.tif
        ├── mosaic_Cellbound2_z5.tif
        ├── mosaic_Cellbound2_z6.tif
        ├── mosaic_Cellbound3_z0.tif
        ├── mosaic_Cellbound3_z1.tif
        ├── mosaic_Cellbound3_z2.tif
        ├── mosaic_Cellbound3_z3.tif
        ├── mosaic_Cellbound3_z4.tif
        ├── mosaic_Cellbound3_z5.tif
        ├── mosaic_Cellbound3_z6.tif
        ├── mosaic_DAPI_z0.tif
        ├── mosaic_DAPI_z1.tif
        ├── mosaic_DAPI_z2.tif
        ├── mosaic_DAPI_z3.tif
        ├── mosaic_DAPI_z4.tif
        ├── mosaic_DAPI_z5.tif
        ├── mosaic_DAPI_z6.tif
        ├── mosaic_PolyT_z0.tif
        ├── mosaic_PolyT_z1.tif
        ├── mosaic_PolyT_z2.tif
        ├── mosaic_PolyT_z3.tif
        ├── mosaic_PolyT_z4.tif
        ├── mosaic_PolyT_z5.tif
        └── mosaic_PolyT_z6.tif

- **Merfish input data directory** → maps to workflow param `vizgen_images`
The user provided input must be provided by the user and wrapped as `LatchDir(latch://...)`.  

#### Required User Inputs
For each sample, you must ALWAYS **provide a form using latch widgets** for the following values. Each one has a clear mapping to a workflow parameter: 

- **algorithm json** → A json file encoding the algorithm parameters `cell_segmentation_algorithm`  (`LatchFile`)
- **Output directory** (where to save Trekker workflow results) → `output_directory` (`LatchOutputDir`)  
  - If not provided, default to:  
    ```python
    LatchDir("latch://38438.account/vizgen_cellsegmentation_outputs")
    ```  
- **Run name** → A string to identify the cell segmentation run and maps to ```run_name```

When you use the w_ldata_picker widget to populate the `output_directory` or `cell_segmentation_algorithm` params, ALWAYS retrieve the LData path string by accessing the widget `.value.path` before passing to LatchFile(...) or LatchDir(...)
                     
Use the code below as a template, that uses w_workflow. The code will generate a “Launch” button. Always make sure you are activating the workflow by ```execution = w.value```, where w is the workflow.

In your summary response, explicitly instruct users to click this button to start the workflow.
                     
#### Example Implementation

```python
from lplots.widgets.workflow import w_workflow
from latch.types import LatchFile, LatchOutputDir

params = {
    "run_name": "run_1",
    "vizgen_images":LatchDir("latch://38438.account/202305010900_U2OS_small_set_VMSC00000"),
    "output_directory": LatchDir("latch://38438.account/vizgen_cellsegmentation_outputs"),
    "cell_segmentation_algorithm": LatchFile("latch://38438.account/202305010900_U2OS_small_set_VMSC00000/algorithm.json")
}

w = w_workflow(
    wf_name="wf.__init__.vizgen_cell_segmentation_wf",
    version="0.0.0-c28480",
    params=params,
    label="Run Cell Segmentation Workflow",
)
execution = w.value
```