cwlVersion: v1.0
class: CommandLineTool
label: Runs fastQC on each fastq file in fastq directory
hints:
  DockerRequirement:
    dockerPull: hubmap/scanpy:latest
baseCommand: /opt/fastqc_wrapper.py

inputs:
  fastq_dir:
    type: Directory
    doc: Directory containing fastq files to be evaluated
    inputBinding:
      position: 1

outputs:
  zipped_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.zip"
    doc: Individual graph files and additional data files containing the raw data from which plots were drawn.

  report_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.html"
    doc: HTML reports with embedded graphs
