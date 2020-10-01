cwlVersion: v1.0
class: CommandLineTool
label: Assay-specific annotation of cell barcodes after quantification
hints:
  DockerRequirement:
    dockerPull: hubmap/scanpy:1.5.1
baseCommand: /opt/annotate_cells.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  h5ad_file:
    type: File
    inputBinding:
      position: 1
  metadata_json:
    type: File?
    inputBinding:
      position: 2
      prefix: '--metadata_json'
outputs:
  annotated_h5ad_file:
    type: File
    outputBinding:
      glob: 'out.h5ad'
