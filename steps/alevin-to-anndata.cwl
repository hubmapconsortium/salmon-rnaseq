cwlVersion: v1.0
class: CommandLineTool
label: Convert Alevin sparse output to anndata.AnnData object, save as h5ad
requirements:
  DockerRequirement:
    dockerPull: hubmap/scrna-analysis:latest
baseCommand: /opt/alevin_to_anndata.py

inputs:
  alevin_dir:
    type: Directory
    inputBinding:
      position: 1
outputs:
  raw_expr_h5ad:
    type: File
    outputBinding:
      glob: raw_expr.h5ad
  expr_h5ad:
    type: File
    outputBinding:
      glob: expr.h5ad
