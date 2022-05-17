cwlVersion: v1.0
class: CommandLineTool
label: RNA velocity analysis via scVelo
hints:
  DockerRequirement:
    dockerPull: hubmap/scrna-analysis:2.1.9b1
baseCommand: /opt/scvelo_analysis.py

inputs:
  spliced_h5ad_file:
    type: File
    inputBinding:
      position: 1
outputs:
  annotated_h5ad_file:
    type: File
    outputBinding:
      glob: scvelo_annotated.h5ad
  embedding_grid_plot:
    type: File
    outputBinding:
      glob: scvelo_embedding_grid.pdf
