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
  fastqc_dir:
    outputSource: fastqc/fastqc_dir
    type: Directory
    label: "Directory of FastQC output files, mirroring input directory structure"
  salmon_output:
      outputSource: salmon/output_dir
      type: Directory
      label: "Full output of `salmon quant`"
  expression_matrix:
    outputSource: make_expression_matrix/expression_matrix
    type: File
    label: "A hd5 file containing transcript by sample matrices of TPM and number of reads"

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
      - quant_files
      - command_info
      - auxiliary_files
    run: steps/bulk-salmon.cwl
    label: "Salmon quant 1.0.0, with index from GRCh38 transcriptome"

  - id: make_expression_matrix
    in:
      - id: quant_files
        source: salmon-bulk/quant_files

    out:
      - expression_matrix
    run: steps/expression-matrix.cwl
    label: "Make expression matrix from quant vectors"
