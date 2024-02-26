cwlVersion: v1.0
class: CommandLineTool
label: Convert Alevin sparse output to anndata.AnnData object, save as h5ad
requirements:
  DockerRequirement:
    dockerPull: hubmap/scrna-analysis:2.1.19
baseCommand: /opt/alevin_to_anndata.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  alevin_dir:
    type: Directory
    inputBinding:
      position: 1
outputs:
  raw_expr_h5ad:
    type: File?
    outputBinding:
      glob: raw_expr.h5ad
  expr_h5ad:
    type: File
    outputBinding:
      glob: expr.h5ad
  genome_build_json:
    type: File
    outputBinding:
      glob: genome_build.json
