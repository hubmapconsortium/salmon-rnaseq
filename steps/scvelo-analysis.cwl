cwlVersion: v1.2
class: CommandLineTool
label: RNA velocity analysis via scVelo
requirements:
  DockerRequirement:
    dockerPull: hubmap/scrna-analysis:2.3.3
  EnvVarRequirement:
    envDef:
      TMPDIR: "/tmp"
baseCommand: /opt/scvelo_analysis.py

inputs:
  spliced_h5ad_file:
    type: File
    inputBinding:
      position: 1

  assay_name:
    type: string
    inputBinding:
      position: 2

outputs:
  annotated_h5ad_file:
    type: File?
    outputBinding:
      glob: scvelo_annotated.h5ad
  embedding_grid_plot:
    type: File?
    outputBinding:
      glob: scvelo_embedding_grid.pdf
