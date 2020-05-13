cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: hubmap/scanpy:20200513-192405

baseCommand: /opt/make_expression_matrix.py
label: Takes gene expression vectors from several bulk RNA samples and makes them into a gene by sample matrix

inputs:
  quant_files:
    type: File[]
    inputBinding:
      position: 0

outputs:
  expression_matrix:
    type: File
    outputBinding:
      glob: expression_matrices.h5
