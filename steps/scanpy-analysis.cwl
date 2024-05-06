cwlVersion: v1.0
class: CommandLineTool
label: Dimensionality reduction and clustering
requirements:
  DockerRequirement:
    dockerPull: hubmap/scrna-analysis:2.2.1
baseCommand: /opt/scanpy_entry_point.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  h5ad_file:
    type: File
    inputBinding:
      position: 1
outputs:
  filtered_data_h5ad:
    type: File
    outputBinding:
      glob: secondary_analysis.h5ad
  dispersion_plot:
    type: File
    outputBinding:
      glob: dispersion_plot.pdf
  umap_plot:
    type: File
    outputBinding:
      glob: umap_by_leiden_cluster.pdf
  spatial_plot:
    type: File?
    outputBinding:
      glob: spatial_pos_by_leiden_cluster.pdf
  umap_density_plot:
    type: File
    outputBinding:
      glob: umap_embedding_density.pdf
  marker_gene_plot_t_test:
    type: File
    outputBinding:
      glob: marker_genes_by_cluster_t_test.pdf
  marker_gene_plot_logreg:
    type: File
    outputBinding:
      glob: marker_genes_by_cluster_logreg.pdf
