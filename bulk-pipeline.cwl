#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: bulk scRNA-seq pipeline using Salmon
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

steps:

  - id: fastqc
    in:
      - id: fastq_dir
        source: fastq_dir
    out:
      - zipped_files
      - report_files
    run: steps/fastqc.cwl
    label: "Run fastqc on all fastq files in fastq directory"

  - id: salmon-bulk
    in:
      - id: fastq_dir
        source: fastq_dir
      - id: threads
        source: threads
    out:
      - quant_dir
      - quant_sfs
      - command_info
      - auxiliary_files
    run: steps/bulk-salmon.cwl
    label: "Salmon quant 1.0.0, with index from GRCh38 transcriptome"

  - id: make_expression_matrix
    in:
      - id: quant_dir
        source: salmon-bulk/quant_dir

    out:
      - expression_matrix
    run: steps/expression-matrix.cwl
    label: "Make expression matrix from quant vectors"
