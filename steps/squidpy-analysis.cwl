cwlVersion: v1.0
class: CommandLineTool
label: Dimensionality reduction and clustering
requirements:
  DockerRequirement:
    dockerPull: hubmap/squidpy-analysis:2.2.1
baseCommand: /opt/squidpy_entry_point.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  h5ad_file:
    type: File
    inputBinding:
      position: 1
  img_dir:
    type: Directory?
    inputBinding:
      position: 2
outputs:
  squidpy_annotated_h5ad:
    type: File?
    outputBinding:
      glob: squidpy_annotated.h5ad
  neighborhood_enrichment_plot:
    type: File?
    outputBinding:
      glob: neighborhood_enrichment.pdf
  co_occurrence_plot:
    type: File?
    outputBinding:
      glob: co_occurrence.pdf
  spatial_plot:
    type: File?
    outputBinding:
      glob: spatial_scatter.pdf
  interaction_matrix_plot:
    type: File?
    outputBinding:
      glob: interaction_matrix.pdf
  centrality_scores_plot:
    type: File?
    outputBinding:
      glob: centrality_scores.pdf
  ripley_plot:
    type: File?
    outputBinding:
      glob: ripley.pdf
