cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: hubmap/salmon-grch38:1.5
baseCommand: /opt/salmon_wrapper.py
label: Run Salmon Alevin tool on FASTQ input

# arguments are hardcoded in salmon_wrapper.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  adj_fastq_dir:
    type: Directory
    inputBinding:
      position: 1
  orig_fastq_dirs:
    type: Directory[]
    inputBinding:
      position: 2
  threads:
    type: int
    inputBinding:
      position: 3
      prefix: "--threads"

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: salmon_out
