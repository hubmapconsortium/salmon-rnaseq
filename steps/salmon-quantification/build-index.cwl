cwlVersion: v1.0
class: CommandLineTool
label: Assay-specific adjustment of cell barcodes
requirements:
  DockerRequirement:
    dockerPull: hubmap/salmon-grch38:visium-ffpe
    baseCommand: /opt/build_salmon_index.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  fastq_dir:
    type: Directory[]
    inputBinding:
      position: 1
outputs:
  salmon_index:
    type: File?
    outputBinding:
      glob: 'visium_index"
