#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: scRNA-seq pipeline using Salmon and Alevin
#requirements:
#  ScatterFeatureRequirement: {}
#  SubworkflowFeatureRequirement: {}
inputs:
  fastq_dir:
    label: "Directory containing FASTQ files"
    type: Directory
  threads:
    label: "Number of threads for Salmon"
    type: int
    default: 1
outputs:
  salmon_output:
    outputSource: salmon/output_dir
    type: Directory
    label: "Full output of `salmon alevin`"
  count_matrix:
    outputSource: alevin_to_anndata/h5ad_file
    type: File
    label: "Unfiltered count matrix from Alevin, converted to H5AD"
  zipped_files:
    outputSource: fastqc/zipped_files
    type:
      type: array
      items: File
    label: "Individual graph files and additional data files containing the raw data from which plots were drawn."
  report_files:
    outputSource: fastqc/report_files
    type:
      type: array
      items: File
    label: "HTML reports with embedded graphs"
  qc_results:
    outputSource: qc_checks/qc_results
    type: File
    label: "Quality control metrics"
  umap_pdf:
    outputSource: dim_reduce_cluster/umap_pdf
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
    outputSource: cluster_diffexpr/marker_gene_plot_t_test
    type: File
    label: "Cluster marker genes, t-test"
  marker_gene_plot_logreg:
    outputSource: cluster_diffexpr/marker_gene_plot_logreg
    type: File
    label: "Cluster marker genes, logreg method"
steps:
  - id: salmon
    in:
      - id: fastq_dir
        source: fastq_dir
      - id: threads
        source: threads
    out:
      - quant_mat
      - quant_mat_cols
      - quant_mat_rows
      - quant_tier_mat
    run: steps/salmon.cwl
    label: "Salmon Alevin, with index from GRCh38 transcriptome"
  - id: fastqc
    in:
      - id: fastq_dir
        source: fastq_dir
    out:
      - zipped_files
      - report_files
    run: steps/fastqc.cwl
    label: "Run fastqc on all fastq files in fastq directory"
  - id: alevin_to_anndata
    in:
      - id: quant_mat
        source: salmon/quant_mat
      - id: quant_mat_cols
        source: salmon/quant_mat_cols
      - id: quant_mat_rows
        source: salmon/quant_mat_rows
    out:
      - h5ad_file
    run: steps/alevin-to-anndata.cwl
    label: "Convert Alevin output to AnnData object in h5ad format"
  - id: scanpy_analysis
    in:
      - id: h5ad_file
        source: alevin_to_anndata/h5ad_file
    out:
      - qc_results
      - filtered_data
      - umap_pdf
      - marker_gene_plot_t_test
      - marker_gene_plot_logreg
    run: steps/scanpy-analysis.cwl
    label: "Secondary analysis via ScanPy"
