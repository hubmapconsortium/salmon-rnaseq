cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: mruffalo/salmon-grch38:20200213
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

  quant_dir:
    type: Directory
    outputBinding:
      glob: out/salmon_quant/

  quant_sfs:
    type: array
    items:
      type: array
      items: File
    outputBinding:
      glob: /out/salmon_quant/*quant.sf

  command_info:
    type: File
    outputBinding:
      glob: /out/salmon_quant/cmd_info.json

  auxiliary_files:
    type: File
    outputBinding:
      glob: /out/salmon_quant/aux_files.tar.gz
