#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.1
label: scRNA-seq pipeline using Salmon and Alevin
#requirements:
#  ScatterFeatureRequirement: {}
#  SubworkflowFeatureRequirement: {}
inputs:
  fastq_r1:
    label: "FASTQ file"
    type: File
  fastq_r2:
    label: "FASTQ file"
    type: File
  dataset_id:
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
  count_matrix:
    outputSource: alevin_to_anndata/h5ad_file
    type: File
    label: "Unfiltered count matrix from Alevin, converted to H5AD"

steps:
  correct_barcodes:
    in:
      fastq_r1:
        source: fastq_r1
      fastq_r2:
        source: fastq_r2
    out: [barcode_umi_fastq, transcript_fastq]
    run: steps/correct-barcodes.cwl
  salmon:
    in:
      barcode_umi_fastq:
        source: correct_barcodes/barcode_umi_fastq
      transcript_fastq:
        source: correct_barcodes/transcript_fastq
      threads:
        source: threads
    out: [output_dir]
    run: steps/salmon.cwl
    label: "Salmon Alevin, with index from GRCh38 transcriptome"
  alevin_to_anndata:
    in:
      alevin_dir:
        source: salmon/output_dir
      dataset_id:
        source: dataset_id
    out: [h5ad_file]
    run: steps/alevin-to-anndata.cwl
    label: "Convert Alevin output to AnnData object in h5ad format"
