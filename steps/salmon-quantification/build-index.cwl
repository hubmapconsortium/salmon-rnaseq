cwlVersion: v1.0
class: CommandLineTool
label: Assay-specific adjustment of cell barcodes
requirements:
  DockerRequirement:
    dockerPull: hubmap/salmon-grch-38:latest
baseCommand: /opt/build_salmon_index.py

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
  salmon_index:
    type: File?
    outputBinding:
      glob: 'adj_fastq'
