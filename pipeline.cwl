#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: scRNA-seq pipeline using Salmon and Alevin
requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
inputs:
  fastq_dir:
    label: "Directory containing FASTQ files"
    type: Directory[]
  assay:
    label: "scRNA-seq assay"
    type: string
  threads:
    label: "Number of threads for Salmon"
    type: int
    default: 1
outputs:
  salmon_output:
    outputSource: salmon/output_dir
    type: Directory
    label: "Full output of `salmon alevin`"
  count_matrix_h5ad:
    outputSource: annotate_cells/annotated_h5ad_file
    type: File
    label: "Unfiltered count matrix from Alevin, converted to H5AD, spliced and unspliced counts"
  count_matrix_zarr:
    outputSource: annotate_cells/annotated_zarr_dir
    type: Directory
    label: "Unfiltered count matrix from Alevin, converted to H5AD, spliced and unspliced counts"
  raw_count_matrix:
    outputSource: alevin_to_anndata/raw_expr_h5ad
    type: File
    label: "Unfiltered count matrix from Alevin, converted to H5AD, with intronic regions"
  fastqc_dir:
    outputSource: fastqc/fastqc_dir
    type: Directory[]
    label: "Directory of FastQC output files, mirroring input directory structure"
  qc_results:
    outputSource: scanpy_analysis/qc_results
    type: File
    label: "Quality control metrics"
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
  filtered_data_h5ad:
    outputSource: scanpy_analysis/filtered_data_h5ad
    type: File
    label: Full data set of filtered results
    doc: >-
      Full data set of filtered results: expression matrix, coordinates in
      dimensionality-reduced space (PCA and UMAP), cluster assignments via
      the Leiden algorithm, and marker genes for one cluster vs. rest
  filtered_data_zarr:
    outputSource: scanpy_analysis/filtered_data_zarr
    type: Directory
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
    type: File
    label: "scVelo-annotated h5ad file, including cell RNA velocity"
  scvelo_annotated_zarr:
    outputSource: scvelo_analysis/annotated_zarr_dir
    type: Directory
    label: "scVelo-annotated zarr directory, including cell RNA velocity"
  scvelo_embedding_grid_plot:
    outputSource: scvelo_analysis/embedding_grid_plot
    type: File
    label: "scVelo velocity embedding grid plot"
  scvelo_embedding_stream_plot:
    outputSource: scvelo_analysis/embedding_stream_plot
    type: File?
    label: "scVelo velocity embedding stream plot"
steps:
  adjust_barcodes:
    in:
      fastq_dir:
        source: fastq_dir
      assay:
        source: assay
    out: [adj_fastq_dir, metadata_json]
    run: steps/adjust-barcodes.cwl
  salmon:
    in:
      orig_fastq_dirs:
        source: fastq_dir
      adj_fastq_dir:
        source: adjust_barcodes/adj_fastq_dir
      assay:
        source: assay
      threads:
        source: threads
    out:
      - output_dir
    run: steps/salmon.cwl
    label: "Salmon Alevin, with index from GRCh38 transcriptome"
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
  alevin_to_anndata:
    in:
      alevin_dir:
        source: salmon/output_dir
    out:
      - expr_h5ad
      - raw_expr_h5ad
    run: steps/alevin-to-anndata.cwl
    label: "Convert Alevin output to AnnData object in h5ad format"
  annotate_cells:
    in:
      assay:
        source: assay
      h5ad_file:
        source: alevin_to_anndata/expr_h5ad
      metadata_json:
        source: adjust_barcodes/metadata_json
    out:
      - annotated_h5ad_file
      - annotated_zarr_dir
    run: steps/annotate-cells.cwl
  scanpy_analysis:
    in:
      h5ad_file:
        source: annotate_cells/annotated_h5ad_file
    out:
      - qc_results
      - filtered_data_h5ad
      - filtered_data_zarr
      - umap_plot
      - marker_gene_plot_t_test
      - marker_gene_plot_logreg
      - dispersion_plot
      - umap_density_plot
    run: steps/scanpy-analysis.cwl
    label: "Secondary analysis via ScanPy"
  scvelo_analysis:
    in:
      spliced_h5ad_file:
        source: annotate_cells/annotated_h5ad_file
    out:
      - annotated_h5ad_file
      - annotated_zarr_dir
      - embedding_grid_plot
      - embedding_stream_plot
    run: steps/scvelo-analysis.cwl
    label: "RNA velocity analysis via scVelo"
