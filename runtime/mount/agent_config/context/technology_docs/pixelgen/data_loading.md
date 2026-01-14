<goal>
Load Pixelator PNA data from workflow outputs (PXL files) and prepare AnnData object with sample and condition annotations.
</goal>

<method>
1/ **Import Required Libraries**
   - Import Pixelator read function:
     ```python
     from pixelator import read_pna as read
     ```

2/ **Identify and Download PXL Files**
   - Identify PXL files from workflow entries (files ending with ".pxl")
   - Download each PXL file locally
   - Extract sample names by removing suffixes in order:
     - ".layout.dataset.pxl"
     - ".layout.pxl"
     - ".dataset.pxl"
     - ".pxl"
   - Store sample names and local file paths

3/ **Load and Aggregate PNA Data**
   - Use read_pna to aggregate PXL files:
     ```python
     pg_data_combined = read(local_paths)
     ```
   - Extract AnnData object:
     ```python
     adata = pg_data_combined.adata()
     ```

4/ **Add Condition Annotations**
   - Option A: Automatic condition assignment based on sample name patterns
     - Allow user to specify pattern matching rules (e.g., samples containing "PHA" → "PHA stimulated", others → "Resting")
     - Apply condition labels:
       ```python
       adata.obs["condition"] = "condition_1"
       adata.obs.loc[adata.obs["sample"].astype(str).str.contains("pattern", na=False), "condition"] = "condition_2"
       ```
   - Option B: Manual condition assignment
     - Allow user to manually assign conditions to each sample
     - Create condition mapping and apply to adata.obs
</method>

<workflows>
</workflows>

<library>
</library>

<self_eval_criteria>
</self_eval_criteria>
