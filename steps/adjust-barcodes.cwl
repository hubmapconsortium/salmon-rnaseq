cwlVersion: v1.0
class: CommandLineTool
label: Assay-specific adjustment of cell barcodes
hints:
  DockerRequirement:
    dockerPull: hubmap/scanpy:1.4
baseCommand: /opt/adjust_barcodes.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  fastq_dir:
    type: Directory[]
    inputBinding:
      position: 1
outputs:
  adj_fastq_dir:
    type: Directory
    outputBinding:
      glob: 'adj_fastq'
  metadata_json:
    type: File?
    outputBinding:
      glob: 'metadata.json'
