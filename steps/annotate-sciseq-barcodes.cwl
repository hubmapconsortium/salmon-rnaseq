cwlVersion: v1.0
class: CommandLineTool
label: Annotate sci-seq barcodes
hints:
  DockerRequirement:
    dockerPull: hubmap/scanpy-sciseq:1.3-sci
baseCommand: /opt/annotate_sciseq_barcodes.py

inputs:
  h5ad_file:
    type: File
    inputBinding:
      position: 1
outputs:
  annotated_h5ad_file:
    type: File
    outputBinding:
      glob: 'out.h5ad'
