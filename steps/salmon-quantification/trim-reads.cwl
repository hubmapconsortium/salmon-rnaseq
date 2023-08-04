cwlVersion: v1.0
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/scrna-trim-reads:2.1.13
baseCommand: /opt/trim_reads.py
label: Trim FASTQ files

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
  trimmed_fastq_dir:
    type: Directory
    outputBinding:
      glob: trimmed
