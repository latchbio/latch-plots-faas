## Vizgen Merfish Analysis Workflow

This is the **authoritative step-by-step pipeline** for **Vizgen Merfish** experiments.  Always follow steps **in order**. 

1. **Experiment Setup** — **ALWAYS** check if the users provide a h5ad file and a latch directory encoding spatial images. 
    - If the user provided only the spatial images ask them to provide the h5ad file.
    - If the user provided only the h5ad file, ASK them for the spatial directory. 
    - Check if the spatial directory directory contains `pmtiles` and if they are raw images, inform the user. Spatial coordinates are typically found in **spatial** or **X_spatial**
2. **Data Loading** — load the `h5ad data` using **Scanpy**.  Ensure the input is a valid **AnnData (H5AD)** object. - Always use `w_h5` that uses anndata and **spatial directory** as inputs
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
11. **Cell Type Annotation** — Prompt the user to input tissue and organism. Then Assign **biological meaning** to clusters using marker-gene dictionaries to identify cell types. Make a dot plot of cell-types to identify top markers using scanpy.
---

✅ **General Rules**  
- Always produce **visual outputs** at every stage. 
- Use `w_plot()` and `w_table()` to display plots and tables on the UI
- Use `w_text_output()` to display progress and information to users
- Make plots with plotly unless explicitly specified
- Apply sensible defaults
- Add well formatted markdown
- Prioritize **transparency**, **reproducibility**, and **interactivity** throughout the analysis.

Some pointers on steps.

### Experiment Setup

- You must ask users to confirm whether their experiment is vizgen Merfish
- If the user provides an `h5ad` file, prompt the user to input a directory containining spatial files

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

- **Merfish input data directory** → maps to workflow param `vizgen_images`.
The user provided input must be provided by the user and wrapped as `LatchDir(latch://...)`.  

### Required User Inputs
For each sample, you must ALWAYS **provide a form using latch widgets** for the following values. Each one has a clear mapping to a workflow parameter: 

#### Task: Create a New algorithm.json File
Overview
    - Always generate a new algorithm.json file based on user input.
    - Ignore any existing algorithm.json files in the directory.

#### Steps:
1. Look for **manifest.json** file in the **Merfish Input Directory**
    - The manifest.json file looks similar to, 
    ```
    { 
        "microns_per_pixel": 0.108, 
        "mosaic_width_pixels": 3,953, 
        "mosaic_height_pixels": 3,960, 
        "bbox_microns": [ 
            -4.274, 
            -6.39, 
            422.5, 
            421.2
        ], 
        "hor_num_tiles_box": 2, 
        "vert_num_tiles_box": 2, 
        "mosaic_files": [ 
            { 
                "stain": "Cellbound3", 
                "z": 0, 
                "file_name": "mosaic_Cellbound3_z0.tif"
            }, 
            { 
                "stain": "Cellbound2", 
                "z": 0, 
                "file_name": "mosaic_Cellbound2_z0.tif"
            },
        ...
        ]
    }
    ```
    - Detect and list all zoom levels (`z`) and stains (`stain`) available in the **manifest.json** file.
    - Display a summary of the detected zoom levels and stains. **Ensure there is atleast one zoom level and one stain**.

2. Number of Segmentation Tasks
    - Use an **lplots text input widget** to ask the user for the number of segmentation tasks.
    - If the user enters N, create N segmentation configuration sections.

3. Configure Each Segmentation Task
    - **USE A FOR LOOP** to match the number of segmentation tasks and define parameters for each segmentation task
    - For each segmentation task 1..N:
        - Display the segmentation task number using **w_text_output**
        - Display available zoom levels using **lplots widgets**(default: all selected).
        - Allow selection of the segmentation algorithm:
            - Options: Cellpose, Cellpose2, Watershed
            - Default: Cellpose2
            - If Cellpose or Cellpose2 is selected:
                - Use **lplots widgets** to assign red, green, and blue stains (from detected stains).
                - Use an **lplots widgets** to choose the segmentation model:
                    - Options: cyto2, nuclei
                    - Default: cyto2
                - This selection configures the "model" field in segmentation_parameters.
            - If  Watershed is selected:
                - Do not include channel map fields.
        - Always set the following fields in segmentation_parameters:
        ```
        {
            "nuclear_channel": "all",
            "entity_fill_channel": "all"
        }
        ```

4. Dynamic Updates
    - If the user changes the number of segmentation tasks, automatically regenerate the segmentation configuration
to match the updated number using **lplots signals**. 

5. Save Algorithm File
    - Use an **lplots widgets** to let the user specify the file name for the new JSON document.
    - Use an **ldatapicker widget** to let the user choose the directory in LData to save the file.
    - Combine both inputs to save the final algorithm.json file to LData.

#### Important Notes
    - Always create a new algorithm.json regardless of existing files.
    - The file content must reflect the latest user inputs.
    - The workflow must handle multiple segmentation configurations dynamically and clearly.
The **algorithm.json** file should look similar to 
    ```
    {
    "experiment_properties": {
        "all_z_indexes": [
        0,
        1,
        2,
        3,
        4,
        5,
        6
        ],
        "z_positions_um": [
        1.5,
        3,
        4.5,
        6,
        7.5,
        9,
        10.5
        ]
    },
    "segmentation_tasks": [
        {
        "task_id": 0,
        "segmentation_family": "Cellpose2",
        "entity_types_detected": [
            "cell"
        ],
        "z_layers": [
            0, 1, 2, 3, 4, 5, 6
        ],
        "segmentation_properties": {
            "model": "cyto2",
            "model_dimensions": "2D",
            "custom_weights": "CP_20230830_093420",
            "channel_map": {
            "red": "Cellbound1",
            "green": "Cellbound3",
            "blue": "DAPI"
            }
        },
        "task_input_data": [
            {
            "image_channel": "Cellbound1",
            "image_preprocessing": [
                {
                "name": "normalize",
                "parameters": {
                    "type": "CLAHE",
                    "clip_limit": 0.01,
                    "filter_size": [
                    100,
                    100
                    ]
                }
                }
            ]
            },
            {
            "image_channel": "Cellbound3",
            "image_preprocessing": [
                {
                "name": "normalize",
                "parameters": {
                    "type": "CLAHE",
                    "clip_limit": 0.01,
                    "filter_size": [
                    100,
                    100
                    ]
                }
                }
            ]
            },
            {
            "image_channel": "DAPI",
            "image_preprocessing": [
                {
                "name": "normalize",
                "parameters": {
                    "type": "CLAHE",
                    "clip_limit": 0.01,
                    "filter_size": [
                    100,
                    100
                    ]
                }
                }
            ]
            }
        ],
        "segmentation_parameters": {
            "nuclear_channel": "DAPI",
            "entity_fill_channel": "all",
            "diameter": 70,
            "flow_threshold": 0.95,
            "cellprob_threshold": -5.5,
            "minimum_mask_size": 500
        },
        "polygon_parameters": {
            "simplification_tol": 2,
            "smoothing_radius": 10,
            "minimum_final_area": 500
        }
        },
        {
        "task_id": 1,
        "segmentation_family": "Cellpose2",
        "entity_types_detected": [
            "cell"
        ],
        "z_layers": [
            3
        ],
        "segmentation_properties": {
            "model": "nuclei",
            "model_dimensions": "2D",
            "custom_weights": null,
            "channel_map": {
            "red": "",
            "green": "",
            "blue": "DAPI"
            }
        },
        "task_input_data": [
            {
            "image_channel": "DAPI",
            "image_preprocessing": [
                {
                "name": "normalize",
                "parameters": {
                    "type": "CLAHE",
                    "clip_limit": 0.01,
                    "filter_size": [
                    100,
                    100
                    ]
                }
                }
            ]
            }
        ],
        "segmentation_parameters": {
            "nuclear_channel": "all",
            "entity_fill_channel": "DAPI",
            "diameter": 55,
            "flow_threshold": 0.8,
            "cellprob_threshold": -3.0,
            "minimum_mask_size": 500
        },
        "polygon_parameters": {
            "simplification_tol": 2,
            "smoothing_radius": 10,
            "minimum_final_area": 500
        }
        }
    ],
    "segmentation_task_fusion": {
        "entity_fusion_strategy": "harmonize",
        "fused_polygon_postprocessing_parameters": {
        "min_distance_between_entities": 1,
        "min_final_area": 500
        }
    },
    "output_files": [
        {
        "entity_types_output": [
            "cell"
        ],
        "files": {
            "run_on_tile_dir": "result_tiles/",
            "mosaic_geometry_file": "mosaic_space.parquet",
            "micron_geometry_file": "micron_space.parquet",
            "cell_metadata_file": "cell_metadata.csv"
        }
        }
    ]
    }
    ```
   
- **algorithm json** → The json file encoding the algorithm parameters produced in the previous step `cell_segmentation_algorithm`  (`LatchFile`)
- **Output directory** (where to save cell segmentation workflow results) → `output_directory` (`LatchOutputDir`)  
  - If not provided, default to:  

    ```python
    LatchDir("latch://38438.account/vizgen_cellsegmentation_outputs")
    ```
  
- **Run name** → A string to identify the cell segmentation run and maps to `run_name`

When you use the `w_ldata_picker` widget to populate the `output_directory` or `cell_segmentation_algorithm` params, ALWAYS retrieve the LData path string by accessing the widget `.value.path` before passing to LatchFile(...) or LatchDir(...)
                     
Use the code below as a template, that uses `w_workflow`. The code will generate a “Launch” button. Always make sure you are activating the workflow by `execution = w.value`, where w is the workflow.

In your summary response, explicitly instruct users to click this button to start the workflow.

#### Input directory to workflow

**ALWAYS USE THE USER PROVIDED DIRECTORY AS THE INPUT to `vizgen_images`**   
The input should point to the top-level dataset directory, not its subfolders.
DO:
    ```
    LatchDir("latch://.../202305010900_U2OS_small_set_VMSC00000")
    ```
DO NOT:
    ```
    LatchDir("latch://.../202305010900_U2OS_small_set_VMSC00000/region_0")
    ```

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
