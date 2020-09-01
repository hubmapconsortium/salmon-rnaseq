cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: hubmap/salmon-snare-grch38:1.3.1-snare
baseCommand: /opt/salmon_wrapper.py
label: Run Salmon Alevin tool on FASTQ input

# arguments are hardcoded in salmon_wrapper.py

inputs:
  threads:
    type: int
    inputBinding:
      position: 0
      prefix: "--threads"
  barcode_umi_fastq:
    type: File
    inputBinding:
      position: 1
  transcript_fastq:
    type: File
    inputBinding:
      position: 2

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: out
