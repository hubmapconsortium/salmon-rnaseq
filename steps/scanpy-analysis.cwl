cwlVersion: v1.0
class: CommandLineTool
label: Dimensionality reduction and clustering
hints:
  DockerRequirement:
    dockerPull: hubmap/scrna-analysis:latest
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
  filtered_data_zarr:
    type: Directory
    outputBinding:
      glob: secondary_analysis.zarr
  dispersion_plot:
    type: File
    outputBinding:
      glob: dispersion_plot.pdf
  umap_plot:
    type: File
    outputBinding:
      glob: umap_by_leiden_cluster.pdf
  slideseq_plot:
    type: File?
    outputBinding:
      glob: spatial_pos_by_leiden_cluster.pdf
  umap_density_plot:
    type: File
    outputBinding:
      glob: umap_embedding_density.pdf
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
