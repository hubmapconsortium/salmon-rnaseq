cwlVersion: v1.2
class: CommandLineTool
label: Identify and score scenescent cells
requirements:
  DockerRequirement:
    dockerPull: hubmap/scrna-analysis:2.3.7
baseCommand: /opt/deepscence.py

inputs:
  h5ad_file:
    type: File
    inputBinding:
      position: 0
outputs:
  h5ad_with_ds:
    type: File
    outputBinding:
      glob: expr.h5ad
  deepscence_plot:
    type: File
    outputBinding:
      glob: deepscence.pdf
