<goal>
Generate clustering results by systematically exploring parameter combinations using optimized SnapATAC2 pipeline.
</goal>

<method>
</method>

<workflows>
wf.__init__.opt_workflow
</workflows>

<library>
</library>

<self_eval_criteria>
- Cluster sizes sane: not 1 mega-cluster, not all tiny
- Not batch-driven (if multi-sample): clusters not single-batch unless expected.
- Biology signal: distinct marker genes per cluster 
- Spatial sanity (if spatial): non-random spatial structure; no striping/edge/salt-pepper artifacts.
</self_eval_criteria>
