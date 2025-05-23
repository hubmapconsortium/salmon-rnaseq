cwlVersion: v1.2
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/salmon-grcm39:latest
  ResourceRequirement:
    ramMin: 28672
baseCommand: /opt/salmon_wrapper.py
label: Run Salmon Alevin tool on FASTQ input

# arguments are hardcoded in salmon_wrapper.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  trimmed_fastq_dir:
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
  expected_cell_count:
    type: int?
    inputBinding:
      position: 4
      prefix: "--expected-cell-count"
  keep_all_barcodes:
    type: boolean?
    inputBinding:
      position: 5
      prefix: "--keep-all-barcodes"
  organism:
    type: string?
    inputBinding:
      position: 6
      prefix: "--organism"

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: salmon_out
