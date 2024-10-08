#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.2
label: scRNA-seq pipeline using Salmon and Alevin
requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
inputs:
  fastq_dir:
    label: "Directory containing FASTQ files"
    type: Directory[]
  img_dir:
    label: "Directory containing TIFF image data (for Visium assay)"
    type: Directory?
  metadata_dir:
    label: "Directory containing metadata, including gpr slide file (for Visium assay)"
    type: Directory?
  assay:
    label: "scRNA-seq assay"
    type: string
  threads:
    label: "Number of threads for Salmon"
    type: int
    default: 1
  expected_cell_count:
    type: int?
  keep_all_barcodes:
    type: boolean?
  organism:
    type: string?
outputs:
  salmon_output:
    outputSource: salmon_quantification/salmon_output
    type: Directory
    label: "Full output of `salmon alevin`"
  count_matrix_h5ad:
    outputSource: salmon_quantification/count_matrix_h5ad
    type: File
    label: "Unfiltered count matrix from Alevin, converted to H5AD, spliced and unspliced counts"
  raw_count_matrix:
    outputSource: salmon_quantification/raw_count_matrix
    type: File?
    label: "Unfiltered count matrix from Alevin, converted to H5AD, with intronic counts as separate columns"
  genome_build_json:
    outputSource: salmon_quantification/genome_build_json
    type: File
    label: "Genome build information in JSON format"
  fastqc_dir:
    outputSource: fastqc/fastqc_dir
    type: Directory[]
    label: "Directory of FastQC output files, mirroring input directory structure"
  scanpy_qc_results:
    outputSource: compute_qc_results/scanpy_qc_results
    type: File
    label: "Quality control metrics from Scanpy"
  qc_report:
    outputSource: compute_qc_results/qc_metrics
    type: File
    label: "Quality control report in JSON format"
  dispersion_plot:
    outputSource: scanpy_analysis/dispersion_plot
    type: File
    label: "Gene expression dispersion plot"
  umap_plot:
    outputSource: scanpy_analysis/umap_plot
    type: File
    label: "UMAP dimensionality reduction plot"
  umap_density_plot:
    outputSource: scanpy_analysis/umap_density_plot
    type: File
    label: "UMAP dimensionality reduction plot, colored by cell density"
  spatial_plot:
    outputSource: scanpy_analysis/spatial_plot
    type: File?
    label: "Slide-seq bead plot, colored by Leiden cluster"
  filtered_data_h5ad:
    outputSource: scanpy_analysis/filtered_data_h5ad
    type: File
    label: Full data set of filtered results
    doc: >-
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
  scvelo_annotated_h5ad:
    outputSource: scvelo_analysis/annotated_h5ad_file
    type: File?
    label: "scVelo-annotated h5ad file, including cell RNA velocity"
  scvelo_embedding_grid_plot:
    outputSource: scvelo_analysis/embedding_grid_plot
    type: File?
    label: "scVelo velocity embedding grid plot"
  squidpy_annotated_h5ad:
    outputSource: squidpy_analysis/squidpy_annotated_h5ad
    type: File?
  neighborhood_enrichment_plot:
    outputSource: squidpy_analysis/neighborhood_enrichment_plot
    type: File?
  co_occurrence_plot:
    outputSource: squidpy_analysis/co_occurrence_plot
    type: File?
  interaction_matrix_plot:
    outputSource: squidpy_analysis/interaction_matrix_plot
    type: File?
  centrality_scores_plot:
    outputSource: squidpy_analysis/centrality_scores_plot
    type: File?
  ripley_plot:
    outputSource: squidpy_analysis/ripley_plot
    type: File?
  squidpy_spatial_plot:
    outputSource: squidpy_analysis/spatial_plot
    type: File?
steps:
  salmon_quantification:
    in:
      fastq_dir:
        source: fastq_dir
      img_dir:
        source: img_dir
      metadata_dir:
        source: metadata_dir
      assay:
        source: assay
      threads:
        source: threads
      expected_cell_count:
        source: expected_cell_count
      keep_all_barcodes:
        source: keep_all_barcodes
      organism:
        source: organism
    out:
      - salmon_output
      - count_matrix_h5ad
      - raw_count_matrix
      - genome_build_json
    run: steps/salmon-quantification.cwl
  fastqc:
    scatter: [fastq_dir]
    scatterMethod: dotproduct
    in:
      fastq_dir:
        source: fastq_dir
      threads:
        source: threads
    out:
      - fastqc_dir
    run: steps/fastqc.cwl
    label: "Run fastqc on all fastq files in fastq directory"
  scanpy_analysis:
    in:
      assay:
        source: assay
      h5ad_file:
        source: salmon_quantification/count_matrix_h5ad
    out:
      - filtered_data_h5ad
      - umap_plot
      - marker_gene_plot_t_test
      - marker_gene_plot_logreg
      - dispersion_plot
      - umap_density_plot
      - spatial_plot
    run: steps/scanpy-analysis.cwl
    label: "Secondary analysis via ScanPy"
  scvelo_analysis:
    in:
      spliced_h5ad_file:
        source: salmon_quantification/count_matrix_h5ad
      assay_name:
        source: assay
    out:
      - annotated_h5ad_file
      - embedding_grid_plot
    run: steps/scvelo-analysis.cwl
    label: "RNA velocity analysis via scVelo"
  squidpy_analysis:
    in:
      assay:
        source: assay
      h5ad_file:
        source: scanpy_analysis/filtered_data_h5ad
      img_dir:
        source: img_dir
    out:
      - squidpy_annotated_h5ad
      - neighborhood_enrichment_plot
      - co_occurrence_plot
      - interaction_matrix_plot
      - ripley_plot
      - centrality_scores_plot
      - spatial_plot
    run: steps/squidpy-analysis.cwl
    label: "Spatial analysis via SquidPy"
  compute_qc_results:
    in:
      assay:
        source: assay
      primary_matrix_path:
        source: salmon_quantification/count_matrix_h5ad
      secondary_matrix_path:
        source: scanpy_analysis/filtered_data_h5ad
      salmon_dir:
        source: salmon_quantification/salmon_output
    out:
      - scanpy_qc_results
      - qc_metrics
    run: steps/compute-qc-metrics.cwl
    label: "Compute QC metrics"