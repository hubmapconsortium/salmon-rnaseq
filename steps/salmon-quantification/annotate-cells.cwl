cwlVersion: v1.0
class: CommandLineTool
label: Assay-specific annotation of cell barcodes after quantification
requirements:
  DockerRequirement:
    dockerPull: hubmap/scrna-analysis:2.2
baseCommand: /opt/annotate_cells.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  h5ad_file:
    type: File
    inputBinding:
      position: 1
  orig_fastq_dirs:
    type: Directory[]
    inputBinding:
      position: 2
  img_dir:
    type: Directory?
    inputBinding:
      position: 3
      prefix: '--img_dir'
  metadata_dir:
    type: Directory?
    inputBinding:
      position: 4
      prefix: '--metadata_dir'
  metadata_json:
    type: File?
    inputBinding:
      position: 5
      prefix: '--metadata_json'
outputs:
  annotated_h5ad_file:
    type: File
    outputBinding:
      glob: 'expr.h5ad'
