cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: mruffalo/salmon-grch38:20200213
baseCommand: /opt/make_expression_matrix.py
label: Takes gene expression vectors from several bulk RNA samples and makes them into a gene by sample matrix

inputs:
  quant_sf_dir:
    type: Directory
    inputBinding:
      position: 1

outputs:
  expression_matrix:
    type: File
    outputBinding:
      glob: out/expression_matrices.h5
