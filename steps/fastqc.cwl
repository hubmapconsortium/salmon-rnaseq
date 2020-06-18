cwlVersion: v1.0
class: CommandLineTool
label: Runs fastQC on each fastq file in fastq directory
hints:
  DockerRequirement:
    dockerPull: hubmap/scanpy:1.2
baseCommand: /opt/fastqc_wrapper.py

inputs:
  fastq_dir:
    type: Directory
    doc: Directory containing fastq files to be evaluated
    inputBinding:
      position: 1
  threads:
    type: int
    doc: The number of threads to use for fastqc
    inputBinding:
      position: 2

outputs:
  fastqc_dir:
    type: Directory
    outputBinding:
      glob: "fastqc_output"
    doc: Individual graph files and additional data files containing the raw data from which plots were drawn.
