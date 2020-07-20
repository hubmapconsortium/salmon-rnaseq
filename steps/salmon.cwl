cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: hubmap/salmon-sciseq-grch38:1.3-sci
baseCommand: /opt/salmon_wrapper.py
label: Run Salmon Alevin tool on FASTQ input

# arguments are hardcoded in salmon_wrapper.py

inputs:
  threads:
    type: int
    inputBinding:
      position: 0
      prefix: "--threads"
  adj_fastq_dir:
    type: Directory
    inputBinding:
      position: 1

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: out
