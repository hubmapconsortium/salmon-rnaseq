cwlVersion: v1.0
class: CommandLineTool
label: Compute QC metrics
requirements:
  DockerRequirement:
    dockerPull: hubmap/scrna-analysis:2.1
baseCommand: /opt/compute_qc_metrics.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  h5ad_primary:
    type: File
    inputBinding:
      position: 1
  h5ad_secondary:
    type: File
    inputBinding:
      position: 2
  salmon_dir:
    type: Directory
    inputBinding:
      position: 3
outputs:
  scanpy_qc_results:
    type: File
    outputBinding:
      glob: qc_results.hdf5
  qc_metrics:
    type: File
    outputBinding:
      glob: qc_results.json
