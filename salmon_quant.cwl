#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.1
label: scRNA-seq pipeline using Salmon and Alevin
#requirements:
#  ScatterFeatureRequirement: {}
#  SubworkflowFeatureRequirement: {}
inputs:
  fastq_dir:
    label: "FASTQ directory"
    type: Directory[]
  threads:
    label: "Number of threads for Salmon"
    type: int
    default: 1
outputs:
  output_dir:
    outputSource: salmon/output_dir
    type: Directory
    label: "Full output of `salmon alevin`"

steps:
  correct_barcodes:
    in:
      fastq_dir:
        source: fastq_dir
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
