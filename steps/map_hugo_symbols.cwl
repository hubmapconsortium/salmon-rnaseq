cwlVersion: v1.2
class: CommandLineTool
label: Identify and score scenescent cells
requirements:
  DockerRequirement:
    dockerPull: hubmap/mygene:latest
baseCommand: /opt/map_hugo_symbols.py

inputs:
  h5ad_file:
    type: File
    inputBinding:
      position: 0
outputs:
  h5ad_with_hugo_symbols:
    type: File
    outputBinding:
      glob: secondary_analysis.h5ad
