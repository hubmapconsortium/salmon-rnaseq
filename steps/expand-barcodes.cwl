cwlVersion: v1.0
class: CommandLineTool
label: Expand sci-seq barcodes
hints:
  DockerRequirement:
    dockerPull: hubmap/scanpy-sciseq:latest
baseCommand: /opt/expand_sciseq_barcodes.py

inputs:
  fastq_dir:
    type: Directory
    inputBinding:
      position: 1
outputs:
  adj_fastq_dir:
    type: Directory
    outputBinding:
      glob: 'adj_fastq'
