cwlVersion: v1.0
class: CommandLineTool
label: Convert Alevin sparse output to anndata.AnnData object, save as h5ad
hints:
  DockerRequirement:
    dockerPull: hubmap/scanpy:20200513-152658
baseCommand: /opt/alevin_to_anndata.py

inputs:
  alevin_dir:
    type: Directory
    inputBinding:
      position: 1
outputs:
  h5ad_file:
    type: File
    outputBinding:
      glob: out.h5ad
