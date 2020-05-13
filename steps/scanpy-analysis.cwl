cwlVersion: v1.0
class: CommandLineTool
label: Dimensionality reduction and clustering
hints:
  DockerRequirement:
    dockerPull: hubmap/scanpy:latest
baseCommand: /opt/scanpy_entry_point.py

arguments:
  - dim_reduce_cluster
inputs:
  h5ad_file:
    type: File
    inputBinding:
      position: 1
outputs:
  filtered_data:
    type: File
    outputBinding:
      glob: cluster_marker_genes.h5ad
  umap_pdf:
    type: File
    outputBinding:
      glob: umap_by_leiden_cluster.pdf
  qc_results:
    type: File
    outputBinding:
      glob: qc_results.hdf5
  marker_gene_plot_t_test:
    type: File
    outputBinding:
      glob: marker_genes_by_cluster_t_test.pdf
  marker_gene_plot_logreg:
    type: File
    outputBinding:
      glob: marker_genes_by_cluster_logreg.pdf
