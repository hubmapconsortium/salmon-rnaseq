cwlVersion: v1.1
class: Workflow

inputs:
  fastq_dir:
    label: "Directory containing FASTQ files"
    type: Directory
  threads:
    label: "Number of threads for Salmon"
    type: int
    default: 1
outputs:
  count_matrix:
    outputSource: concatenate_h5ad/h5ad_file
    type: File
    label: "Unfiltered count matrix from Alevin, converted to H5AD"
  qc_results:
    outputSource: scanpy_analysis/qc_results
    type: File
    label: "Quality control metrics"
  umap_pdf:
    outputSource: scanpy_analysis/umap_pdf
    type: File
    label: "UMAP dimensionality reduction plot"
  filtered_data:
    outputSource: scanpy_analysis/filtered_data
    type: File
    label: >-
      Full data set of filtered results: expression matrix, coordinates in
      dimensionality-reduced space (PCA and UMAP), cluster assignments via
      the Leiden algorithm, and marker genes for one cluster vs. rest
  marker_gene_plot_t_test:
    outputSource: scanpy_analysis/marker_gene_plot_t_test
    type: File
    label: "Cluster marker genes, t-test"
  marker_gene_plot_logreg:
    outputSource: scanpy_analysis/marker_gene_plot_logreg
    type: File
    label: "Cluster marker genes, logreg method"

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

steps:
  gather_sequence_bundles:
    run: gather_sequence_bundles.cwl
    in:
      sequence_directory: fastq_dir
    out:
      [fastq_r1, fastq_r2, dataset_ids]

  salmon_quant:
    scatter: [fastq_r1, fastq_r2, dataset_id]
    scatterMethod: dotproduct
    run: salmon_quant.cwl
    in:
      fastq_r1: gather_sequence_bundles/fastq_r1
      fastq_r2: gather_sequence_bundles/fastq_r2
      dataset_id: gather_sequence_bundles/dataset_ids
      threads: threads
    out: [count_matrix, salmon_output]

  concatenate_h5ad:
    run: steps/concatenate-h5ad.cwl
    in:
      h5ad_files: salmon_quant/count_matrix
    out: [h5ad_file]

  scanpy_analysis:
    in:
      h5ad_file:
        source: concatenate_h5ad/h5ad_file
    out:
      - qc_results
      - filtered_data
      - umap_pdf
      - marker_gene_plot_t_test
      - marker_gene_plot_logreg
    run: steps/scanpy-analysis.cwl
    label: "Secondary analysis via ScanPy"
