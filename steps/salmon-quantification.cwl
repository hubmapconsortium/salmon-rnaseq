#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Salmon quantification, FASTQ -> H5AD count matrix
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
  expected_cell_count:
    type: int?
  keep_all_barcodes:
    type: boolean?
outputs:
  salmon_output:
    outputSource: salmon/output_dir
    type: Directory
    label: "Full output of `salmon alevin`"
  count_matrix_h5ad:
    outputSource: annotate_cells/annotated_h5ad_file
    type: File
    label: "Unfiltered count matrix from Alevin, converted to H5AD, spliced and unspliced counts"
  raw_count_matrix:
    outputSource: alevin_to_anndata/raw_expr_h5ad
    type: File
    label: "Unfiltered count matrix from Alevin, converted to H5AD, with intronic counts as separate columns"
  genome_build_json:
    outputSource: alevin_to_anndata/genome_build_json
    type: File
    label: "Genome build information in JSON format"

steps:
  adjust_barcodes:
    in:
      fastq_dir:
        source: fastq_dir
      assay:
        source: assay
    out: [adj_fastq_dir, metadata_json]
    run: salmon-quantification/adjust-barcodes.cwl
  trim_reads:
    in:
      orig_fastq_dirs:
        source: fastq_dir
      adj_fastq_dir:
        source: adjust_barcodes/adj_fastq_dir
      assay:
        source: assay
      threads:
        source: threads
    out: [trimmed_fastq_dir]
    run: salmon-quantification/trim-reads.cwl
  salmon:
    in:
      orig_fastq_dirs:
        source: fastq_dir
      trimmed_fastq_dir:
        source: trim_reads/trimmed_fastq_dir
      assay:
        source: assay
      threads:
        source: threads
      expected_cell_count:
        source: expected_cell_count
      keep_all_barcodes:
        source: keep_all_barcodes
    out:
      - output_dir
    run: salmon-quantification/salmon.cwl
    label: "Salmon Alevin, with index from GRCh38 transcriptome"
  alevin_to_anndata:
    in:
      alevin_dir:
        source: salmon/output_dir
    out:
      - expr_h5ad
      - raw_expr_h5ad
      - genome_build_json
    run: salmon-quantification/alevin-to-anndata.cwl
    label: "Convert Alevin output to AnnData object in h5ad format"
  annotate_cells:
    in:
      orig_fastq_dirs:
        source: fastq_dir
      assay:
        source: assay
      h5ad_file:
        source: alevin_to_anndata/expr_h5ad
      metadata_json:
        source: adjust_barcodes/metadata_json
    out:
      - annotated_h5ad_file
    run: salmon-quantification/annotate-cells.cwl
