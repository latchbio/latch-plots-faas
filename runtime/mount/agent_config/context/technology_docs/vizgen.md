# Vizgen Merfish Analysis Workflow
This document provides **step-by-step pipeline** for **Vizgen Merfish** experiments.  
Always follow steps **in order**. ALWAYS use **lplots widgets such as `w_text_input`,`w_text_output`, `w_checkbox`, `w_select`, `w_multi_select`, `w_radio_button_group` to configure parameters for various analysis steps. 

---

## Mandatory Package Redemption

Invoke the `redeem_package` tool to install required Vizgen tools into the workspace.
  - `package_code`: `2b6fe6b03ceef66f2474c604b6cb3a2c5858ce8029f9480bdafd5ebeb60d08a8`
  - `package_version_id`: `400`

## **DATA LOADING**

Follow instructions at `vizgen/vizgen_data_loading.md`

---

## **QUALITTY CONTROL**

Follow instructions at `vizgen/vizgen_QC.md`

## **PREPROCESSING**

**ALWAYS** Check if your adata has > 100,000 cells,
IF TRUE
## **ALWAYS FOLLOW INSTRUCTIONS FROM `rapids.md`**
Procced to **Spatial Analysis Step** after the rapids workflow finishes.

ELSE Skip steps 4-10:
## **ALWAYS FOLLOW INSTRUCTIONS FROM `vizgen_scanpy_preprocessing.md`**
Procced to **Spatial Analysis Step** after this step

---

## **SPATIAL ANALYSIS (using Squidpy)**

Follow instructions at `vizgen/vizgen_spatial_analysis.md`

---

## **SECONDARY ANALYSIS**

Follow instructions at `vizgen/vizgen_secondary_analysis.md`

---

## Launch Cell segmentation Workflow

Follow instructions at `vizgen/vizgen_cell_segmentation.md`

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



