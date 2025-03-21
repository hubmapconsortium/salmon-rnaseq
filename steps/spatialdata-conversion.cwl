cwlVersion: v1.2
class: CommandLineTool
label: Dimensionality reduction and clustering
requirements:
  DockerRequirement:
    dockerPull: hubmap/scrna-spatialdata-conversion:latest
baseCommand: /opt/squidpy_entry_point.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  h5ad_file:
    type: File
    inputBinding:
      position: 1
  img_dir:
    type: Directory?
    inputBinding:
      position: 2
outputs:
  sdata_zarr:
    type: Directory?
    outputBinding:
      glob: squidpy_annotated.h5ad
