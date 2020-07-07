cwlVersion: v1.0
class: CommandLineTool
label: Concatenate H5AD files
hints:
  DockerRequirement:
    dockerPull: mruffalo/scanpy-snare:latest
baseCommand: /opt/concatenate_h5ad.py

inputs:
  h5ad_files:
    type: File[]
    inputBinding:
      position: 1
outputs:
  h5ad_file:
    type: File
    outputBinding:
      glob: out.h5ad
