#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: scRNA-seq pipeline using Salmon and Alevin
requirements:
  SubworkflowFeatureRequirement: {}
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
  fastqc_dir:
    outputSource: fastqc/fastqc_dir
    type: Directory
    label: "Directory of FastQC output files, mirroring input directory structure"
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
steps:
  salmon:
    in:
      fastq_dir:
        source: fastq_dir
      threads:
        source: threads
    out:
      - output_dir
    run: salmon_quant.cwl
    label: "Salmon Alevin, with index from GRCh38 transcriptome"
  fastqc:
    in:
      fastq_dir:
        source: fastq_dir
      threads:
        source: threads
    out:
      - fastqc_dir
    run: steps/fastqc.cwl
    label: "Run fastqc on all fastq files in fastq directory"
  alevin_to_anndata:
    in:
      alevin_dir:
        source: salmon/output_dir
    out:
      - h5ad_file
    run: steps/alevin-to-anndata.cwl
    label: "Convert Alevin output to AnnData object in h5ad format"
  scanpy_analysis:
    in:
      h5ad_file:
        source: alevin_to_anndata/h5ad_file
    out:
      - qc_results
      - filtered_data
      - umap_pdf
      - marker_gene_plot_t_test
      - marker_gene_plot_logreg
    run: steps/scanpy-analysis.cwl
    label: "Secondary analysis via ScanPy"
