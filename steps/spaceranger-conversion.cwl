cwlVersion: v1.0
class: CommandLineTool
label: Dimensionality reduction and clustering
requirements:
  DockerRequirement:
    dockerPull: hubmap/squidpy-analysis:2.1.16
baseCommand: /opt/spaceranger_conversion.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  spaceranger_dir:
    type: Directory
    inputBinding:
      position: 1
outputs:
  raw_spaceranger_h5ad:
    type: File?
    outputBinding:
      glob: raw_spaceranger.h5ad
  filtered_spaceranger_h5ad:
    type: File?
    outputBinding:
      glob: filtered_spaceranger.h5ad
