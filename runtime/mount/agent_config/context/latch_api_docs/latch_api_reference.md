# Latch API Reference

## Widgets

**Core Concepts**

- Every widget is a Python object. Assign to a variable; the live user selection is at `.value`.
- Widgets render in **call order**.
- **Reactive execution:** any cell that touches `some_widget.value` re-runs when that value changes. To read once without reactivity, use `some_widget._signal.sample()`.

---

### w_ldata_picker

**Import:** `from lplots.widgets.ldata import w_ldata_picker`

**When to use:** Select files or directories from Latch Data

**Arguments:**
- `label` (str, required): Label for the picker
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `appearance` (FormInputAppearance, optional): Styling options
- `default` (str, optional): Default latch:// path
- `required` (bool, optional): Whether selection is required. Default: False
- `file_type` (Literal["file", "dir", "any"], optional): Filter by type. Default: "any"
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.ldata import w_ldata_picker
from latch.ldata.path import LPath
import pandas as pd

csv = w_ldata_picker(
    label="CSV",
    default="latch:///path/file.csv",
    file_type="file"
)

if csv.value is not None:
    lp: LPath = csv.value
    fname = lp.name()
    suffix = Path(fname).suffix if fname else ""
    local_p = Path(f"{lp.node_id()}{suffix}")
    lp.download(local_p, cache=True)
    df = pd.read_csv(local_p)

```

---

### w_ldata_browser

**Import:** `from lplots.widgets.ldata import w_ldata_browser`

**When to use:** Display contents of a specific Latch Data directory for browsing

**Arguments:**
- `label` (str, required): Label for the browser
- `dir` (str | LPath, required): Latch Data directory path to browse
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `appearance` (FormInputAppearance, optional): Styling options
- `required` (bool, optional): Whether selection is required. Default: False
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.ldata import w_ldata_browser
from latch.ldata.path import LPath

browser = w_ldata_browser(label="Browse Files", dir="latch:///my_dir")
# browser.value returns the LPath to the directory
```

---

### w_datasource_picker

**Import:** `from lplots.widgets.datasource import w_datasource_picker`

**When to use:** Select from multiple data source types (Latch Data, Registry, DataFrame, Viewer)

**Arguments:**
- `label` (str, required): Label for the picker
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `default` (DataSourceValue, optional): Default datasource selection (discriminated union)
- `appearance` (FormInputAppearance, optional): Styling options
- `required` (bool, optional): Whether selection is required. Default: False
- `key` (str, optional): Unique widget identifier

**Default value types:**
- `{"type":"ldata","node_id":str}` - Latch Data node
- `{"type":"registry","table_id":str}` - Registry table
- `{"type":"dataframe","key":str}` - Kernel DataFrame
- `{"type":"viewer","viewer_id":str}` - Viewer data

**Example:**
```python
from lplots.widgets.datasource import w_datasource_picker

ds = w_datasource_picker(
    label="Datasource",
    default={"type":"ldata","node_id":"95902"}
)
# ds.value returns a pandas DataFrame
```

---

### w_registry_table_picker

**Import:** `from lplots.widgets.registry import w_registry_table_picker`

**When to use:** Select a Registry table to load as DataFrame

**Arguments:**
- `label` (str, required): Label for the picker
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `default` (str, optional): Default table ID
- `appearance` (FormInputAppearance, optional): Styling options
- `required` (bool, optional): Whether selection is required. Default: False
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.registry import w_registry_table_picker
from latch.registry.table import Table

t = w_registry_table_picker(label="Registry table")
if (tid := t.value):
    df = Table(id=tid).get_dataframe()
```

---

### w_registry_table

**Import:** `from lplots.widgets.registry import w_registry_table`

**When to use:** Display a specific Registry table with row selection capability

**Arguments:**
- `label` (str, required): Label for the table viewer
- `table_id` (str, required): Registry table ID to display
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `default` (str, optional): Default value
- `appearance` (FormInputAppearance, optional): Styling options
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.registry import w_registry_table

rt = w_registry_table(label="Sample Data", table_id="12345")
# rt.value returns {"table": Table, "selected_rows": list[Record]}
```

---

### w_dataframe_picker

**Import:** `from lplots.widgets.dataframe import w_dataframe_picker`

**When to use:** Select from DataFrames currently in the notebook kernel

**Arguments:**
- `label` (str, required): Label for the picker
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `appearance` (FormInputAppearance, optional): Styling options
- `required` (bool, optional): Whether selection is required. Default: False
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.dataframe import w_dataframe_picker

df_picker = w_dataframe_picker(label="Select DataFrame")
selected_df = df_picker.value  # Returns the selected DataFrame
```

---

### w_text_input

**Import:** `from lplots.widgets.text import w_text_input`

**When to use:** Get text input from user

**Arguments:**
- `label` (str, required): Label for the input
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `default` (str, optional): Default text value
- `appearance` (FormInputAppearance, optional): Styling options (placeholder, detail, help_text, error_text, description)
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.text import w_text_input

name = w_text_input(label="Name", default="Alice")
user_name = name.value  # Returns the entered text
```

---

### w_text_output

**Import:** `from lplots.widgets.text import w_text_output`

**When to use:** Display text output to user with optional message box styling

**Arguments:**
- `content` (str, required): Text content to display
- `appearance` (TextOutputAppearance, optional): Styling with message_box option
- `key` (str, optional): Unique widget identifier

**Message box types:** "danger", "info", "success", "warning", "primary", "neutral"

**Example:**
```python
from lplots.widgets.text import w_text_output

w_text_output(
    content="Hi " + name.value,
    appearance={"message_box": "success"}
)
```

---

### w_select

**Import:** `from lplots.widgets.select import w_select`

**When to use:** Single selection from a list of options

**Arguments:**
- `label` (str, required): Label for the selector
- `options` (Iterable[str | int | float | bool | datetime], required): Available options
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `default` (str | int | float | bool | datetime, optional): Default selected value
- `appearance` (FormInputAppearance, optional): Styling options
- `required` (bool, optional): Whether selection is required. Default: False
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.select import w_select

sel = w_select(
    label="Choose",
    options=["a", "b", "c"],
    default="a",
    required=False
)
choice = sel.value  # Returns selected option
```

---

### w_multi_select

**Import:** `from lplots.widgets.multiselect import w_multi_select`

**When to use:** Multiple selection from a list of options

**Arguments:**
- `label` (str, required): Label for the selector
- `options` (Iterable[str | int | float | bool | datetime], required): Available options
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `default` (Iterable[str | int | float | bool | datetime], optional): Default selected values
- `appearance` (FormInputAppearance, optional): Styling options
- `required` (bool, optional): Whether selection is required. Default: False
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.multiselect import w_multi_select

ms = w_multi_select(
    label="Tags",
    options=["alpha", "bravo", "charlie"],
    default=["alpha"],
    required=False
)
selected = ms.value  # Returns list of selected options
```

---

### w_radio_group

**Import:** `from lplots.widgets.radio import w_radio_group`

**When to use:** Single selection using radio buttons

**Arguments:**
- `label` (str, required): Label for the radio group
- `options` (Iterable[str | int | float | bool | datetime], required): Available options
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `default` (str | int | float | bool | datetime, optional): Default selected value
- `appearance` (FormInputAppearance, optional): Styling options
- `required` (bool, optional): Whether selection is required. Default: False
- `direction` (Literal["horizontal", "vertical"], optional): Layout direction. Default: "horizontal"
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.radio import w_radio_group

rg = w_radio_group(
    label="One",
    options=[1, 2, 3],
    default=1,
    direction="horizontal",
    required=False
)
selected = rg.value  # Returns selected option
```

---

### w_checkbox

**Import:** `from lplots.widgets.checkbox import w_checkbox`

**When to use:** Boolean checkbox input

**Arguments:**
- `label` (str, required): Label for the checkbox
- `default` (bool, optional): Default checked state. Default: False
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `appearance` (CheckboxInputAppearance, optional): Styling with error_text and description fields
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.checkbox import w_checkbox

cb = w_checkbox(
    label="Flag",
    default=False,
    appearance={"error_text": "Required", "description": "Check to enable"}
)
is_checked = cb.value  # Returns bool
```

---

### w_number_slider_input

**Import:** `from lplots.widgets.slider import w_number_slider_input`

**When to use:** Numeric input with slider interface

**Arguments:**
- `label` (str, required): Label for the slider
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `default` (int | float, optional): Default value
- `min` (int | float, optional): Minimum value
- `max` (int | float, optional): Maximum value
- `step` (int | float, optional): Step increment
- `tooltip_formatter` (Literal["number", "percentage", "currency", "decimal", "integer", "scientific", "bytes"], optional): Format for tooltip
- `scale_type` (Literal["linear", "logarithmic"], optional): Scale type
- `marks` (dict[int | float, str], optional): Custom marks on slider
- `appearance` (FormInputAppearance, optional): Styling options
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.slider import w_number_slider_input

slider = w_number_slider_input(
    label="Threshold",
    default=0.5,
    min=0.0,
    max=1.0,
    step=0.1,
    scale_type="linear",
    tooltip_formatter="percentage"
)
value = slider.value  # Returns numeric value
```

---

### w_range_slider_input

**Import:** `from lplots.widgets.slider import w_range_slider_input`

**When to use:** Select a numeric range with two-handle slider

**Arguments:**
- `label` (str, required): Label for the slider
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `default` (tuple[int, int] | tuple[float, float], optional): Default range values
- `min` (int | float, optional): Minimum value
- `max` (int | float, optional): Maximum value
- `step` (int | float, optional): Step increment
- `tooltip_formatter` (Literal["number", "percentage", "currency", "decimal", "integer", "scientific", "bytes"], optional): Format for tooltip
- `scale_type` (Literal["linear", "logarithmic"], optional): Scale type
- `marks` (dict[int | float, str], optional): Custom marks on slider
- `appearance` (FormInputAppearance, optional): Styling options
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.slider import w_range_slider_input

range_slider = w_range_slider_input(
    label="Range",
    default=(0.2, 0.8),
    min=0.0,
    max=1.0,
    step=0.1
)
range_values = range_slider.value  # Returns tuple (min, max)
```

---

### w_button

**Import:** `from lplots.widgets.button import w_button`

**When to use:** Trigger actions with a button click

**Arguments:**
- `label` (str, required): Button label
- `default` (ButtonWidgetSignalValue, optional): Default button state (internal use)
- `readonly` (bool, optional): Whether button is disabled. Default: False
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.button import w_button

button = w_button(label="Click Button to Run")

if button.value:
    # logic within `if` clause is only run when button is clicked
    print("Button was clicked!")
```

---

### w_plot

**Import:** `from lplots.widgets.plot import w_plot`

**When to use:** Display matplotlib, seaborn, or plotly figures

**Arguments:**
- `label` (str, optional): Label for the plot
- `source` (Figure | SubFigure | Axes | BaseFigure | FacetGrid | PairGrid | JointGrid, optional): Plot object to display
- `appearance` (OutputAppearance, optional): Styling options
- `key` (str, optional): Unique widget identifier

**Note:**
- The plot source must be a named **global variable** for tracking.
- **CRITICAL**: Each `w_plot` must reference **its own unique variable** — reusing a variable will cause all plots to render the same content and create a negative user experience.
- **DO NOT** use `globals()` or dynamic variable naming in loops (e.g., `globals()[f'fig_{i}']`). This does NOT create proper unique variables.
- **CORRECT APPROACH**: Explicitly declare each variable with a unique name (e.g., `fig_plot1`, `fig_plot2`, `fig_plot3`) outside of any loop structure.
- When creating multiple plots, write out each plot creation separately with its own explicit variable name, or build all plots into a list/dict and then create separate variables from that collection.
- For `scanpy` dot or violin plots, you **must first convert the returned object** into **a Matplotlib figure** before passing it to `source` of `w_plot`.

#### Examples

**Basic Usage:**
```python
from lplots.widgets.plot import w_plot
import plotly.express as px

fig_scatter = px.scatter(df, x='x', y='y')
w_plot(label="First Plot", source=fig_scatter)

# Use unique variable for another visualization
fig_bar = px.bar(df, x='x', y='y')
w_plot(label="Second Plot", source=fig_bar)
```

#### ❌ WRONG - Do not do this:
```python
# Using globals() in a loop - plots will NOT render correctly
for i, data in enumerate(datasets):
    fig = px.scatter(data, x='x', y='y')
    globals()[f'fig_{i}'] = fig  # This does NOT work!
    w_plot(label=f"Plot {i}", source=globals()[f'fig_{i}'])
```

#### ✅ RIGHT - Do this instead:

```python
# Explicitly declare each variable
fig_dataset1 = px.scatter(datasets[0], x='x', y='y')
w_plot(label="Dataset 1", source=fig_dataset1)

fig_dataset2 = px.scatter(datasets[1], x='x', y='y')
w_plot(label="Dataset 2", source=fig_dataset2)

fig_dataset3 = px.scatter(datasets[2], x='x', y='y')
w_plot(label="Dataset 3", source=fig_dataset3)
```

```python
# Explicitly convert Scanpy plot objects to matplotlib figures
## Dot plot
dp = sc.pl.dotplot(
    adata,
    var_names=plot_genes,
    groupby='cell_type',
    dendrogram=True,
    return_fig=True,
    show=False,
    figsize=(14, 6),
)
dp.show()
fig_dotplot = dp.fig
w_plot(label="Cell Type Marker Genes Dot Plot", source=fig_dotplot)
```

---

### w_table

**Import:** `from lplots.widgets.table import w_table`

**When to use:** Display pandas DataFrames

**Arguments:**
- `label` (str, optional): Label for the table
- `source` (DataFrame, optional): DataFrame to display
- `appearance` (OutputAppearance, optional): Styling options
- `key` (str, optional): Unique widget identifier

**Note:** The source must be a named global variable for tracking.

**Example:**
```python
from lplots.widgets.table import w_table
import pandas as pd

df = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
w_table(label="Data", source=df)
```

---

### w_h5

**Import:** `from lplots.widgets.h5 import w_h5`

**When to use:** Interactive AnnData/H5AD viewer for single-cell and spatial data

**Arguments:**
- `label` (str, optional): Label for the viewer
- `ann_data` (AnnData, optional): AnnData object to visualize
- `spatial_dir` (LPath, optional): Path to spatial data directory
- `ann_tiles` (LPath, optional): Path to annotation tiles
- `readonly` (bool, optional): Whether viewer is read-only. Default: False
- `appearance` (OutputAppearance, optional): Styling options
- `viewer_presets` (ViewerPreset, optional): Configuration for viewer defaults
- `key` (str, optional): Unique widget identifier

**Viewer presets options:**
- `genes_of_interest` (list[str]): Genes to highlight
- `default_color_by` (ColorByObs | ColorByVar): Default coloring scheme
- `default_obsm_key` (str): Default embedding to display (e.g., "X_umap", "spatial")
- `cell_markers` (CellMarkers): Marker size and opacity settings
- `categorical_color_palette` (list[str]): Colors for categorical data
- `continuous_color_palette` (list[str]): Colors for continuous data

**Example:**
```python
from lplots.widgets.h5 import w_h5
from latch.ldata.path import LPath

viewer = w_h5(
    ann_data=adata,
    readonly=False,
    viewer_presets={
        "genes_of_interest": ["CD3D", "CD4"],
        "default_color_by": {"type": "obs", "key": "cell_type"},
        "default_obsm_key": "X_umap",
        "cell_markers": {"default_size": 3, "default_opacity": 0.8},
        "categorical_color_palette": ["red", "blue"],
        "continuous_color_palette": ["blue", "white", "red"]
    }
)

v = viewer.value
# v["lasso_points"]: list[list[(x,y)]]
# v["lasso_points_obsm"]: str | None
# v["image_alignment_step"]: current step in image alignment workflow
```

---

### w_ann_data

**Import:** `from lplots.widgets.ann_data import w_ann_data`

**When to use:** Display or interact with AnnData objects (lighter alternative to w_h5)

**Arguments:**
- `ann_data` (AnnData, optional): AnnData object to interact with
- `readonly` (bool, optional): Whether widget is read-only. Default: False
- `appearance` (OutputAppearance, optional): Styling options
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.ann_data import w_ann_data

ann_widget = w_ann_data(ann_data=adata, readonly=False)
# ann_widget.value returns the AnnData object
```

---

### w_igv

**Import:** `from lplots.widgets.igv import w_igv, IGVOptions`

**When to use:** Display genomics data (BAM, VCF, BED, BigWig) in IGV browser

**Arguments:**
- `label` (str, optional): Label for the viewer
- `options` (IGVOptions, required): IGV.js configuration dictionary
- `key` (str, optional): Unique widget identifier

**IGVOptions common fields:**
- `genome` (str): Reference genome ID (e.g., "hg38", "mm10")
- `locus` (str | list[str]): Initial genomic locus (e.g., "chr1:100000-200000")
- `tracks` (list[TrackOptions]): List of track configurations
- `showNavigation` (bool): Toggle navigation bar. Default: True
- `showIdeogram` (bool): Toggle ideogram display. Default: True
- `showRuler` (bool): Show base-pair ruler. Default: True
- `readOnly` (bool): Disable editing. Default: False

**Track common fields:**
- `name` (str): Display name of the track
- `type` (str): Track type ("alignment", "variant", "annotation", "wig")
- `url` (str): Path or URL to data file (latch:// or public)
- `indexURL` (str): Path or URL to index file (.bai, .tbi, etc.)
- `color` (str): Track color (name, hex, or rgb)
- `height` (int): Track height in pixels
- `displayMode` (str): "EXPANDED" or "SQUISHED"
- `autoscale` (bool): Automatically adjust y-axis for coverage tracks
- `visibilityWindow` (int): Max visible region in base pairs. Default: 100000

**Example:**
```python
from lplots.widgets.igv import w_igv, IGVOptions
from latch.account import Account

workspace_id = Account.current().id

latch_path = f"latch://{workspace_id}.account/Covid/covid.bam"
index_path = f"latch://{workspace_id}.account/Covid/covid.bam.bai"

options: IGVOptions = {
    "genome": "hg38",
    "locus": "chr1:155,100,000-155,200,000",
    "tracks": [
        {
            "name": "Alignment",
            "type": "alignment",
            "url": latch_path,
            "indexURL": index_path,
            "color": "steelblue",
            "height": 150
        }
    ]
}

w_igv(options=options)
```

---

### w_logs_display

**Import:** `from lplots.widgets.logs import w_logs_display`

**When to use:** Display logs and progress messages

**Arguments:**
- `label` (str, optional): Label for the logs display
- `appearance` (FormInputAppearance, optional): Styling options
- `key` (str, optional): Unique widget identifier

**Note:** Must call `submit_widget_state()` after creating to render properly.

**Example:**
```python
from lplots.widgets.logs import w_logs_display
from lplots import submit_widget_state

w_logs_display()
submit_widget_state()
```

---

### w_workflow

**Import:** `from lplots.widgets.workflow import w_workflow`

**When to use:** Launch a Latch Workflow directly from Plots

**Arguments:**
- `label` (str, required): Button label
- `wf_name` (str, required): Name of the workflow to execute
- `params` (dict, required): Dictionary of input parameters
- `automatic` (bool, required): Launch workflow automatically. Should always be True.
- `key` (str, required): Unique widget identifier -- determines when the workflow will relaunch (subsequent cell runs with the same key will not relaunch the workflow)
- `version` (str, optional): Workflow version; defaults to latest
- `readonly` (bool, optional): Disable button if True. Default: False

**Returns:** `Execution` object or None

**Example:**
```python
from lplots.widgets.workflow import w_workflow

workflow = w_workflow(
    label="Run Analysis",
    automatic=True,
    key="my_analysis_workflow_run_1",
    wf_name="my_analysis_workflow",
    params={
        "input_file": "latch://workspace/data/sample.fastq",
        "output_dir": "latch://workspace/results/"
    },
    version="v1.0"
)

execution = workflow.value

if execution is not None:
  res = await execution.wait()

  if res is not None and res.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
      # inspect workflow outputs for downstream analysis
      workflow_outputs = list(res.output.values())
```

---

### w_row

**Import:** `from lplots.widgets.row import w_row`

**When to use:** Arrange widgets horizontally

**Arguments:**
- `items` (list[BaseWidget], required): List of widgets to arrange
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.row import w_row

w_row(items=[widget1, widget2, widget3])
```

---

### w_column

**Import:** `from lplots.widgets.column import w_column`

**When to use:** Arrange widgets vertically

**Arguments:**
- `items` (list[BaseWidget], required): List of widgets to arrange
- `key` (str, optional): Unique widget identifier

**Example:**
```python
from lplots.widgets.column import w_column

w_column(items=[widget1, widget2, widget3])
```

---

### w_grid

**Import:** `from lplots.widgets.grid import w_grid`

**When to use:** Arrange widgets in a grid with custom spans

**Arguments:**
- `columns` (int, required): Number of columns in the grid
- `rows` (int, optional): Number of rows (auto if not specified)
- `key` (str, optional): Unique widget identifier

**Methods:**
- `add(item, col_span=1, row_span=1)`: Add widget to grid with span

**Example:**
```python
from lplots.widgets.grid import w_grid

with w_grid(columns=12) as g:
    g.add(item=widget1, col_span=4, row_span=1)
    g.add(item=widget2, col_span=8, row_span=1)
```

---

## LPath

### LPath

**Import:** `from latch.ldata.path import LPath`

**When to use:** All remote `latch://` file and directory operations

**Core Rules:**
- Remote (`latch://`) paths → **LPath** only. Local paths → **pathlib.Path** only
- If a widget already returns an **LPath**, use it directly (don't wrap again)
- **Always check for `None`** before using widget values
- **Always explicitly define local destination and cache file downloads**
- **Never** pass LPath directly to libraries expecting local paths; download first

**Valid Path Forms:**
- With domain: `latch://{domain}/path` (two slashes) - e.g., `latch://12345.account/Data/file.csv`
- Without domain: `latch:///path` (three slashes) - relative paths only

**Common Methods:**

**Metadata (lazy loading):**
- `node_id(load_if_missing=True) -> Optional[str]` - Get node ID
- `name(load_if_missing=True) -> Optional[str]` - Get file/directory name
- `type(load_if_missing=True) -> Optional[LDataNodeType]` - Get node type
- `size_recursive(load_if_missing=True) -> Optional[int]` - Get total size including children
- `size(load_if_missing=True) -> Optional[int]` - Get size of this node only
- `content_type(load_if_missing=True) -> Optional[str]` - Get MIME type
- `is_dir(load_if_missing=True) -> bool` - Check if directory
- `fetch_metadata() -> None` - Force network fetch and cache metadata
- `exists(load_if_missing=True) -> bool` - Check if path exists

**Directory Operations:**
- `iterdir() -> Iterator[LPath]` - Iterate over directory contents (non-recursive, network call)
- `mkdirp() -> None` - Create directory and parents if needed
- `rmr() -> None` - Recursively delete (raises if path missing, network call)

**Data Movement:**
- `copy_to(dst: LPath) -> None` - Copy remote to remote
- `upload_from(src: pathlib.Path, show_progress_bar: bool=False) -> None` - Upload local to remote
- `download(dst: Optional[pathlib.Path]=None, show_progress_bar: bool=False, cache: bool=False) -> pathlib.Path` - Download to local

**Path Operations:**
- `LPath / "child.ext"` - Join paths with `/` operator

**DO Examples:**
```python
from latch.ldata.path import LPath
from pathlib import Path
import pandas as pd

# Construct remote path
lp = LPath("latch://XXXXX.account/Data/project/file.csv")

# Join paths
child = LPath("latch://XXXXX.account/Data/project") / "file.csv"

# Get basename (method, not attribute!)
fname = lp.name()
suffix = Path(fname).suffix if fname else ""
local_p = Path(f"{lp.node_id()}{suffix}")

# Download before using with local libraries
# Always cache downloads
lp.download(local_p, cache=True)
df = pd.read_csv(local_p)

# Check if directory
if lp.is_dir():
    for child in lp.iterdir():
        print(child.name())
```

**DON'T Examples:**
```python
# ❌ Using Path for remote
from pathlib import Path
bad = Path("latch://...")

# ❌ Re-wrapping an LPath
bad = LPath(existing_lpath)

# ❌ Calling str() on an LPath
bad = str(lp)

# ❌ os.path.join or f-strings for remote paths
bad = os.path.join("latch://...", "file.csv")
bad = f"{lp}/file.csv"

# ❌ Accessing .name as attribute
bad = lp.name  # use lp.name() instead
```

**Common Workflow Pattern:**
```python
from lplots.widgets.ldata import w_ldata_picker
from latch.ldata.path import LPath
from pathlib import Path
import pandas as pd

# Get file from widget
pick = w_ldata_picker(label="Select file/folder")
if pick.value is None:
    exit(0)

lp: LPath = pick.value  # Already an LPath, use directly

# Get file info
fname = lp.name()
suffix = Path(fname).suffix if fname else ""
local_p = Path(f"{lp.node_id()}{suffix}")

# Always cache downloads
lp.download(local_p, cache=True)

# Use with local libraries
if fname and fname.endswith(".csv"):
    df = pd.read_csv(local_p)
```

**Quick Diagnostics:**
- Getting errors with pandas/seaborn on `latch://`? → Download to local Path first
- Seeing `ValueError` from `LPath(...)`? → You passed an LPath into LPath(); use it directly
- Path fails with domain and `latch:///`? → Use two slashes with domain: `latch://{domain}/path`

---

## Reactivity

### Signal

**Import:** `from lplots.reactive import Signal`

**When to use:** Cross-cell dependencies and reactive execution

**Core Concepts:**
- **Signal**: Reactive value holder
- **Computation node**: Implicit wrapper around each cell that subscribes to signals it reads
- When a signal's value changes, all current subscribers re-run automatically

**Core Rules:**
1. **NEVER** update and call a Signal within the same cell (causes infinite loop)

**API Methods:**

**Create Signal:**
```python
from lplots.reactive import Signal

value = Signal(10)
```

**Read and subscribe:**
```python
x = value()  # subscribes this cell to changes
```

**Read without subscribing:**
```python
x_now = value.sample()  # snapshot; no subscription
```

**Update Signal:**
```python
value(20)  # schedules re-runs of subscribers (applied after cell completes)
```

**Widget Values as Signals:**

Calling `.value` on a widget retrieves the user input and subscribes the current cell:

```python
# Cell 1 – Create widget
text_input = w_text_input(label="Input")

# Cell 2 – Access value (subscribes to Cell 1)
text_input_value = text_input.value
```

**Immutability Rule:**

Signals store a reference; internal mutation is invisible. Copy when updating:

```python
obj = Signal({})

# ❌ No update propagates
obj()["k"] = "v"

# ✅ Correct way
obj({**obj(), "k": "v"})
```

**Avoid Infinite Loops:**

```python
# ❌ Unconditional self-update → loop
x = value()
value(x + 10)

# ✅ Conditional update → loop terminates
x = value()
if x < 50:
    value(x + 10)
```

**Transactions (Update Timing):**

- All `Signal(...)` updates inside a running cell are queued and applied after the cell finishes
- Multiple cells re-running due to one change are grouped into sequential transactions
- Use `sample()` inside the updating cell if you need the previous value during the same run

**Subscriptions and Conditionals:**

A cell's subscription set is rebuilt on every run—it only subscribes to signals it actually read:

```python
if value() > 30:        # subscribes to `value`
    s = value_str()     # subscribes to `value_str` only when branch executes
```

Force subscription without using the value:
```python
value_str()  # unconditional subscribe
```

**Redefining Signals:**

Global variable reassignment updates value and keeps subscribers:
```python
value = Signal(20)  # equivalent to value(20) if value already exists
```

To drop subscribers and recreate:
```python
del value
value = Signal(0)  # fresh signal, no subscribers
```

For local variables, guard creation:
```python
if "value" not in locals():
    value = Signal(10)
```

**Quick Cheatsheet:**

1. `s = Signal(init)` – create
2. `s()` – read and subscribe
3. `s.sample()` – read without subscribing
4. `s(new)` – update (applies after cell completes)
5. Treat values as immutable; write copies, not in-place mutations
6. No unconditional self-updates in subscribed cells; gate updates
7. Subscriptions = signals read on the last run
8. Reactive transactions hide inconsistent intermediate states
9. Global reassignment preserves subscribers; `del` to reset

**Common Patterns:**

**✅ DO: Separate signal roles**

Don't use the same signal for user input and computation output:

```python
# Good:
clusters_request = Signal(None)  # button/UI
clusters_result = Signal(None)   # computed output
# Viz reads only clusters_result
```

**❌ DON'T: Infinite loops**

Never subscribe and update signals in the same cell:

```python
# Single-cell loop
x = s()       # subscribes
s(x + 1)      # write → rerun → loop

# Cross-cell loop
# Cell A
a = A()       # subscribes
B("trigger")

# Cell B
if B():       # subscribes
    A(a + 1)  # write → Cell A runs → Cell B runs → loop
```

**Common pitfall with AnnData:**
```python
# ❌ This WILL cause an infinite loop
adata_signal()
# ... Do scanpy operations ...
adata_signal(adata)
```

**✅ DO: Break cycles**

Use any of these approaches:
- Add a condition/version check
- Use `sample()` instead of `()` when reading before writing
- Write to a different signal

**✅ DO: Name signals by purpose**

Good names prevent accidental reuse and loops:

**Good:**
- `qc_request`, `qc_result`
- `clusters_request`, `clusters_result`
- `adata_raw`, `adata_qc`, `adata_clustered`, `adata_annotated`

**Bad:**
- `tmp`, `data`, `state`

**Best Practice:** Create a new signal for each stage of analysis (e.g., `adata_qc` → `adata_clustered` → `adata_annotated`) instead of constantly overwriting the same signal.
