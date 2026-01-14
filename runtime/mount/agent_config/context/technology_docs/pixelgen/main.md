<pre_analysis_questions>

- Are PXL files available from workflow outputs, or is data already loaded?
- Are there multiple samples/conditions in the dataset?
- Has normalization (DSB/CLR) been applied to the data?
- Are cell types already annotated, or do you need to perform cell type identification?

</pre_analysis_questions>

<plan>
1. Data Loading -> `pixelgen/data_loading.md`
2. Data Integration (QC and Integration Analysis) -> `pixelgen/Data_Integration.md`
3. Normalization (DSB and CLR) -> `pixelgen/normalization.md`
4. Dimensionality Reduction (PCA, UMAP, Leiden Clustering) -> `pixelgen/dimensionality_reduction.md`
5. Cell Type Annotation (Differential Expression + Marker Gating) -> `pixelgen/celltype_annotation.md`
6. Differential Abundance Analysis -> `pixelgen/differential_abundance.md`
</plan>

<self_eval_criteria>
</self_eval_criteria>
