cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: seandonahue5311/bulk:v1.0
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

  quant_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "out/*.sf"

  command_info:
    type:
      type: array
      items: File
    outputBinding:
      glob: "out/*.json"

  auxiliary_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "out/*.gz"
