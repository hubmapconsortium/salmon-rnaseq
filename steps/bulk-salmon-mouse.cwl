cwlVersion: v1.2
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/salmon-grcm39:latest
baseCommand: /opt/bulk_salmon_wrapper.py
label: Run Salmon quant tool on FASTQ input

inputs:
  threads:
    type: int
    inputBinding:
      position: 0
      prefix: "--threads"
  fastq_dir:
    type: Directory
    inputBinding:
      position: 1
  organism:
    type: string?
    inputBinding:
      position: 2
      prefix: "--organism"

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: out
