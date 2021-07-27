cwlVersion: v1.0
class: CommandLineTool
label: Assay-specific adjustment of cell barcodes
requirements:
  DockerRequirement:
    dockerPull: hubmap/scrna-barcode-adj:2.0.9
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
