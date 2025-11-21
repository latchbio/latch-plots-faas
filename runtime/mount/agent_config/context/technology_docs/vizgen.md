# Vizgen Merfish Analysis Workflow
This is the **authoritative step-by-step pipeline** for **Vizgen Merfish** experiments.  
Always follow steps **in order**. ALWAYS use **lplots widgets such as `w_text_input`,`w_text_output`, `w_checkbox`, `w_select`, `w_multi_select`, `w_radio_button_group` to configure parameters for various analysis steps. 

---

## **1. Experiment Setup**
- **ALWAYS** check if users provide both an **H5AD file** and a **Latch directory** containing spatial images.  
  - If the user provided **only spatial images**, **ask** them to provide the H5AD file.  
  - If the user provided **only the H5AD file**, **ask** them for the spatial directory.  
- Check if the spatial directory contains **pmtiles**; if the images are raw, inform the user.  
- Spatial coordinates are typically found in `spatial` or `X_spatial`.

---

## **2. Data Loading**
- Load the H5AD data using **Scanpy**.  
- Ensure the input is a valid **AnnData (H5AD)** object.  
- Always use **`w_h5`** that takes both **AnnData** and **spatial directory** as inputs.

---

## **3. Quality Control & Filtering**

- Use **Scanpy** to compute QC metrics.  
- **ALWAYS** include standard metrics such as:
  - Total counts per cell  
  - Number of detected genes per cell  
  - Mitochondrial gene percentage  
  - **Cell volume distribution** (derived from polygon area or segmentation mask; visualize histogram). You can never ignore this plot. 
- **ALWAYS** make **histograms** of QC metrics.  
- **ALWAYS** expose filtering parameters using `w_text_input`.  
- **ALWAYS** confirm with the user before applying filters using **lplots widgets**.  
- **ALWAYS** report how many cells will be removed before filtering.  
- If the user declines, **revert to the original AnnData**.

---

## **PREPROCESSING**

### **4. Normalization**
- Provide an option between **log1p** and **total count scaling** using **lplots widgets**.  
- Apply **log1p transformation** by default on the QC-filtered dataset.

### **5. Feature Selection**
- Identify **highly variable genes (HVGs)**.  
- Allow the user to choose the number of HVGs using **lplots widgets** (default: **3,000 genes**).  
- Visualize this step with **Scanpy plots**.

### **6–8. Dimensionality Reduction**
- Check if the AnnData object already has **PCA** and **UMAP**.  
- If they exist, **do not recompute** unless the user explicitly chooses to.  
- If recomputation is needed:
  - Run **PCA** first, then **UMAP**.  
  - Expose parameters using **lplots widgets** for:
    - Number of PCs (**default: 10**)  
    - Number of neighbors (**default: 40**)  
- Always visualize both **PCA** and **UMAP** embeddings.

### **9. Clustering**
- Compute neighbors if they don’t exist.  
- Apply **Leiden clustering** on the neighborhood graph.  
- Use **default resolution = 0.3**.  
- Display clusters on both **UMAP** and **spatial embedding**.

---

## **SPATIAL ANALYSIS (using Squidpy)**

All spatial analyses must use **Squidpy** functions (`sq.gr`, `sq.pl`) and rely on spatial coordinates within the AnnData object.

### **10. Centrality Analysis**
- Compute **spatial centrality metrics** (e.g., degree, betweenness, closeness) using `sq.gr.centrality_scores`. If these centrality scores do not exist, compute these sentrality scores and plot them.
- Plot **centrality plots** overlayed on spatial coordinates.  
- Allow the user to select which metric(s) to visualize. 
```python
sq.gr.centrality_scores(adata, cluster_key="leiden")
```

### **11. Co-occurrence Analysis**
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
### **12. Neighborhood Enrichment**
- Compute **neighborhood enrichment** between clusters using `sq.gr.nhood_enrichment`.  
- Visualize using **heatmaps** and **spatial overlays** (`sq.pl.nhood_enrichment`).  
- Optionally display **z-score–normalized enrichment**.

### **13. Ripley’s Statistics**
- Compute **Ripley’s L** statistics using `sq.gr.ripley`.  
- Plot **observed vs. expected K(r)** curves for major clusters.  
- Highlight clusters showing **significant spatial aggregation**.
- Plot the Ripley's L statistic for **EVERY** cluster using plotly. 
```python
mode = "L"
sq.gr.ripley(adata, cluster_key="leiden", mode=mode)

# Access the ripley results - key format is 'leiden_ripley_L'
ripley_key = f"leiden_ripley_{mode}"
if ripley_key not in adata.uns:
    # Try alternative key format
    ripley_key = "leiden_ripley"

ripley_results = adata.uns[ripley_key]
n_clusters = adata.obs['leiden'].nunique()

# Create line plot for all clusters
df = ripley_results["L_stat"]   # tidy dataframe with: bins, leiden, stats

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
    yaxis_title="Ripley's L(r)",
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

### **14. Spatial Autocorrelation (Moran’s I)**
- Compute **Moran’s I** for gene expression spatial autocorrelation using `sq.gr.spatial_autocorr`.  
- Rank genes by Moran’s I score.  
- Display a **table of top 10 spatially autocorrelated genes**.  
- Optionally, show **spatial expression plots** for these genes with `sq.pl.spatial_scatter`.

---

## **SECONDARY ANALYSIS**

### **15. Differential Gene Expression (DGE)**
- Identify **marker genes per cluster** using rank-based DGE tests.  
- Use **t-test** by default; allow user to select other methods using **lplots widgets**.  
- Report **top marker genes per cluster**.  
- Visualize using **Scanpy dot plots**.

### **16. Cell Type Annotation**
- Prompt the user for **tissue** and **organism**.  
- Assign biological meaning to clusters using **marker-gene dictionaries**.  
- Display **dot plots** of cell types and top markers.

## **17. Spatial Domain Detection and Niche Identification

- Rename `adata.obsm["Spatial"]` to `adata.obsm["X_spatial"]`, then save the processed object to `.h5ad`.
    ```python
    if 'Spatial' in adata.obsm.keys():
        adata.obsm['X_spatial'] = adata.obsm['Spatial']

    # Save the processed H5AD file
    output_h5ad_path = "/tmp/merfish_processed_for_domain_detection.h5ad"
    adata.write_h5ad(output_h5ad_path)
    output_lpath = LPath("latch:///xenium/analysis_output/merfish_processed_for_domain_detection.h5ad")
    output_lpath.upload_from(output_h5ad_path)
    local_path = Path(output_h5ad_path)
    output_lpath.upload_from(local_path)
    ```
- Allow the user to a filename and a location to save the input file for the workflow
- Use this file as input for the domain-detection workflow.
- Allow the user to select lambda, runname and output_dir using a form created with **lplots widgets**
```python
params = {
    "input_file": LatchFile("latch:///xenium_processed_for_domain_detection.h5ad"),
    "run_name": "my_run",
    "output_dir": LatchDir("latch:///Domain_detection_output"),
    "lambda_list": "0.5",
}
w = w_workflow(
    wf_name="wf.__init__.domain_detection_wf",
    key="domain_detection_workflow_run_1",
    version=None,
    automatic=True,
    label="Launch Domain Detection Workflow",
    params=params,
)
execution = w.value

if execution is not None:
  res = await execution.wait()

  if res is not None and res.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
      # inspect workflow outputs for downstream analysis
      workflow_outputs = list(res.output.values())
```
---

✅ **General Rules**
- **NEVER** use `w_number_slider_input` widget at any cost
- **DO NOT** delete cells unless explicitly prompted by the user.
- **DO NOT** create duplicated cells.
- **ALWAYS** store each plot in its own variable. For example the qc plots should be in `fig_qc`. Avoid generic names like `fig` to avoid overwriting figures across cells.
- Always use markdown to summarize text output ```w_text_output``` after each analysis step. **DO NOT** use ```print()```.
- Always produce **visual outputs** at every stage. Use ```w_plot()``` and ```w_table()``` to display plots and tables on the UI
- Always use ```w_h5``` for spatial embedding and UMAP. For all other plots, use Plotly instead. Do not use both w_h5 and Plotly for the same visualization to avoid redundancy.
- Always pass a DataFrame object directly to the source argument of ```w_table```. Do not pass a method call or expression (e.g., df.head(), df.round(3), df.sort_values(...)) directly into ```w_table```.
- Apply sensible defaults
- Add well formatted markdown
- Prioritize **transparency**, **reproducibility**, and **interactivity** throughout the analysis.

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
                     

Use the code below as a template, that uses w_workflow. Always use the `automatic` argument or the workflow will not launch. The workflow will launch automatically when the cell is run. Subsequent cell runs with the same key will not relaunch the workflow, so change the key to a new value if you need to relaunch the workflow.
Finally, you need to make sure to wait for the workflow to complete before proceeding. This is included in the code below.

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
    key="vizgen_cellsegmentation_workflow_run_1",
    version="0.0.0-c28480",
    params=params,
    automatic=True,
    label="Run Cell Segmentation Workflow",
)
execution = w.value

if execution is not None:
  res = await execution.wait()

  if res is not None and res.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
      # inspect workflow outputs for downstream analysis
      workflow_outputs = list(res.output.values())
```
