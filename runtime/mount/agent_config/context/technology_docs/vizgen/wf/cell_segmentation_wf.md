<goal>
Run Vizgen MERFISH cell segmentation from raw images.
</goal>

<parameters>
- **Merfish input data directory** → `vizgen_images` (`LatchDir`)
  - User must provide a directory on Latch Data that corresponds to the **top-level Vizgen dataset folder** (containing `region_0`, `region_1`, etc.).
  - The directory typically looks like:
    - `202305010900_U2OS_small_set_VMSC00000/`
      - `region_0/`
        - `202305010900_U2OS_small_set_VMSC00000_region_0.vzg`
        - `detected_transcripts.csv`
        - `images/`
          - `manifest.json`
          - `micron_to_mosaic_pixel_transform.csv`
          - `mosaic_Cellbound*_z*.tif`, `mosaic_DAPI_z*.tif`, `mosaic_PolyT_z*.tif`, etc.
  - **ALWAYS** pass the **top-level dataset directory** to `vizgen_images`, not subfolders:
    - ✅ DO:
      - `LatchDir("latch://.../202305010900_U2OS_small_set_VMSC00000")`
    - ❌ DO NOT:
      - `LatchDir("latch://.../202305010900_U2OS_small_set_VMSC00000/region_0")`
  - This directory must be wrapped as `LatchDir("latch://...")` before being inserted into the `params` dictionary.

- **Task: Create a new `algorithm.json` file** → `cell_segmentation_algorithm` (`LatchFile`)
  - Always create a **new** `algorithm.json` based on user input; **ignore any existing** `algorithm.json` in the dataset.
  - For **each sample**, provide a form using latch widgets to collect all required values. You must **parse user answers**, normalize them, and construct the final JSON and `params` dictionary.

  - **Step 1 – Parse `manifest.json`**
    - Locate `manifest.json` in `vizgen_images/region_*/images`.
    - Parse and extract:
      - `microns_per_pixel`
      - `mosaic_width_pixels`, `mosaic_height_pixels`
      - `bbox_microns`
      - `mosaic_files` entries with fields:
        - `"stain"` (e.g., `"Cellbound1"`, `"Cellbound2"`, `"Cellbound3"`, `"DAPI"`, `"PolyT"`)
        - `"z"` (zoom / z-index)
        - `"file_name"` (e.g., `"mosaic_Cellbound3_z0.tif"`)
    - Detect all available:
      - **z indexes** → `all_z_indexes`
      - **stains** → list of unique `"stain"` values.
    - Display a summary of detected z-levels and stains using lplots widgets.
    - **Ensure there is at least one z-level and at least one stain**; otherwise, stop and report an error.

  - **Step 2 – Ask for number of segmentation tasks**
    - Use an **lplots text input widget** to ask:
      - “How many segmentation tasks do you want to configure?”
    - Let user enter `N`; this defines how many segmentation configuration blocks will be created in `"segmentation_tasks"`.

  - **Step 3 – Configure each segmentation task (1..N)**
    - Use a **for loop** over `task_id` in `range(N)` to build each task.
    - For each task:
      - Display the current task number using `w_text_output`.
      - Show all detected z-levels and allow the user to select which `z_layers` to use (default: **all** z-levels).
      - Allow the user to choose **segmentation algorithm family**:
        - Options: `"Cellpose"`, `"Cellpose2"`, `"Watershed"`
        - Default: `"Cellpose2"`
      - If `"Cellpose"` or `"Cellpose2"` is selected:
        - Use lplots widgets to allow the user to map stains to RGB:
          - `red`, `green`, `blue` ← choose from detected stains.
        - Use an lplots widget to select **segmentation model**:
          - Options: `"cyto2"`, `"nuclei"`
          - Default: `"cyto2"`
        - These values populate:
          - `"segmentation_family"` (e.g., `"Cellpose2"`)
          - `"segmentation_properties.model"` (e.g., `"cyto2"`)
          - `"segmentation_properties.channel_map"` with `"red"`, `"green"`, `"blue"` keys.
      - If `"Watershed"` is selected:
        - Do **not** include a `channel_map` block (no RGB channel mapping).
      - **Always** set the following fields in `segmentation_parameters` for every task:
        ```json
        {
          "nuclear_channel": "all",
          "entity_fill_channel": "all"
        }
        ```
      - Additional numeric parameters (e.g., `diameter`, `flow_threshold`, `cellprob_threshold`, `minimum_mask_size`, polygon settings) may be exposed as widgets as needed and inserted into:
        - `"segmentation_parameters"` and `"polygon_parameters"`.

  - **Step 4 – Dynamic updates to number of tasks**
    - If the user changes the number of segmentation tasks (`N`), use **lplots signals** to:
      - Automatically regenerate the `segmentation_tasks` configuration list so its length matches the updated `N`.
      - Rebuild widgets and internal structures accordingly.

  - **Step 5 – Save `algorithm.json` to LData**
    - Use an lplots widget to capture the **filename** of the algorithm JSON (e.g., `"algorithm_20230501.json"`).
    - Use an **`w_ldata_picker` widget** to let the user choose the output directory in LData.
    - When using `w_ldata_picker`:
      - **Always retrieve the LData path string** via:
        - `picker.value.path`
      - Then combine:
        - `output_path = picker.value.path`
        - `algorithm_path = f"{output_path}/{filename}"`
    - Write the final JSON content to `algorithm_path`.
    - Wrap the resulting path as:
      - `cell_segmentation_algorithm = LatchFile(algorithm_path)`

  - **Final algorithm.json structure**
    - The resulting file must follow the general structure:
      - `"experiment_properties"` with `all_z_indexes` and `z_positions_um`
      - `"segmentation_tasks"` (list of N task blocks)
      - `"segmentation_task_fusion"` (e.g., `entity_fusion_strategy: "harmonize"`)
      - `"output_files"` with fields such as:
        - `"run_on_tile_dir"`, `"mosaic_geometry_file"`, `"micron_geometry_file"`, `"cell_metadata_file"`
    - The exact content reflects **only the latest user inputs** (widgets + manifest-derived information).

- **Output directory** → `output_directory` (`LatchOutputDir` / `LatchDir`)
  - Directory on Latch Data where the cell segmentation workflow outputs will be saved.
  - If the user does not provide one, default to:
    ```python
    LatchDir("latch://38438.account/vizgen_cellsegmentation_outputs")
    ```
  - When using `w_ldata_picker` to populate this parameter:
    - Always pass:
      - `LatchDir(picker.value.path)`  
      not the widget object directly.

- **Run name** → `run_name` (`str`)
  - A short string to identify the cell segmentation run (e.g., `"run_1"`, `"U2OS_region0_seg_v1"`).
  - This is used to distinguish multiple workflow runs in logs and downstream analysis.

- **Workflow launch via `w_workflow`**
  - Use `w_workflow` from `lplots.widgets.workflow` to launch:
    - `wf_name="wf.__init__.vizgen_cell_segmentation_wf"`
  - You must:
    - Set `automatic=True` so the workflow launches when the cell runs.
    - Provide a unique `key` for each new run. Reusing the same `key` will **not** relaunch the workflow.
    - Wait for completion using `await execution.wait()` before proceeding to dependent analysis.
  - The final `params` dictionary must use:
    - `vizgen_images` (top-level dataset directory)
    - `output_directory`
    - `cell_segmentation_algorithm`
    - `run_name`
  exactly as specified above.
</parameters>

<outputs>
</outputs>

<example>
```python
from lplots.widgets.workflow import w_workflow
from latch.types import LatchFile, LatchDir

params = {
    "run_name": "run_1",
    "vizgen_images": LatchDir("latch://38438.account/202305010900_U2OS_small_set_VMSC00000"),
    "output_directory": LatchDir("latch://38438.account/vizgen_cellsegmentation_outputs"),
    "cell_segmentation_algorithm": LatchFile(
        "latch://38438.account/202305010900_U2OS_small_set_VMSC00000/algorithm.json"
    ),
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
</example>
