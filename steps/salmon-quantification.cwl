#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.2
label: Salmon quantification, FASTQ -> H5AD count matrix
requirements:
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}
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
  img_dir:
    type: Directory?
  metadata_dir:
    type: Directory?
  organism:
    type: string?
outputs:
  salmon_output:
    outputSource: [salmon/output_dir, salmon-mouse/output_dir]
    pickValue: first_non_null
    type: Directory
    label: "Full output of `salmon alevin`"
  count_matrix_h5ad:
    outputSource: annotate_cells/annotated_h5ad_file
    type: File
    label: "Unfiltered count matrix from Alevin, converted to H5AD, spliced and unspliced counts"
  raw_count_matrix:
    outputSource: alevin_to_anndata/raw_expr_h5ad
    type: File?
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
      organism:
        source: organism
    out:
      - output_dir
    run: salmon-quantification/salmon.cwl
    when: $(inputs.organism == 'human')
    label: "Salmon Alevin, with index from GRCh38 transcriptome"
  salmon-mouse:
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
      organism:
        source: organism
    out:
      - output_dir
    run: salmon-quantification/salmon-mouse.cwl
    when: $(inputs.organism == 'mouse')
    label: "Salmon Alevin, with index from GRCm39 transcriptome"
  alevin_to_anndata:
    in:
      assay:
        source: assay
      alevin_dir:
        source: [salmon/output_dir, salmon-mouse/output_dir]
        pickValue: first_non_null
      organism:
        source: organism
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
      img_dir:
        source: img_dir
      metadata_dir:
        source: metadata_dir
      metadata_json:
        source: adjust_barcodes/metadata_json
    out:
      - annotated_h5ad_file
    run: salmon-quantification/annotate-cells.cwl
