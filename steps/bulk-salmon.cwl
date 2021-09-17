cwlVersion: v1.0
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/salmon-grch38:2.1.1
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

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: out
