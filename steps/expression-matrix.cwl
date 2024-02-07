cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: hubmap/scrna-analysis:2.1.18

baseCommand: /opt/make_expression_matrix.py
label: Takes gene expression vectors from several bulk RNA samples and makes them into a gene by sample matrix

inputs:
  quant_dir:
    type: Directory
    inputBinding:
      position: 0

outputs:
  expression_matrix:
    type: File
    outputBinding:
      glob: expression_matrices.h5
