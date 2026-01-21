<goal>
Identify cell types.
</goal>

<method>

**Load priors**: Load `technology_docs/cell_types.json` and select matching entry by (organism, tissue, context). Perform annotation using ONLY `expected_major_types` in priors.

**Trekker (single-cell)**: 
- Use cluster-based cell type annotation. Extract markers per cluster and assign labels based on marker enrichment and biological context.

**Seeker (beads = 1–2 cells per spot)**:
- **Reference-free**: Annotate each **spot** based on marker gene expression.
- **Reference-based**: Use deconvolution via `technology_docs/tools/cell2location.md`. If a single-cell reference is required, use `technology_docs/tools/curate_sc_reference.md`.
</method>

<workflows>
</workflows>

<library>
</library>

<self_eval_criteria>
- Follow `technology_docs/evals/cell_typing.md`
</self_eval_criteria>
