version 1.0
## Use double '#' for workflow-level comments
## This workflow implements a one-task workflow

import "steps/fastqc.wdl" as FastQC
import "steps/salmon-quantification.wdl" as SalmonQuantification

workflow RunSalmonRNAseq {
    meta {
        description: "The SalmonRNAseq pipeline processes"
        author: "Xiang Li"
    }

	input {
        Int threads = 1
        Int mem_gb
        File fastq1
        File fastq2
        String species
        String assay
        String protocol
        String? run_id
        Int? expected_cell_count
        File? limits
        Boolean? keep_all_barcodes
        String? email
	}

    parameter_meta {
        fastq1: "fastq1"
        fastq2: "fastq2"
        threads: "# threads for compute"
        mem_gb: "# memory for compute"
        expected_cell_count: "expected_cell_count"
        keep_all_barcodes: "keep_all_barcodes"
        email: "email"
        species: "species"
        protocol: "protocol"
        run_id: "run_id"
        assay: "assay"
        limits: "limits"
    }

    scatter (fastq in [fastq1, fastq2]) {
        call FastQC.fastqc as fastqc {
            input:
                fastqs = [fastq],
                threads = threads,
                mem_gb = mem_gb,
                limits = limits
        }
    }

    call SalmonQuantification.salmon_quantification as SalmonQuantificationCall {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            threads = threads,
            mem_gb = mem_gb,
            assay = assay,
            species = species,
            protocol = protocol,
            expected_cell_count = expected_cell_count,
            keep_all_barcodes = keep_all_barcodes
    }

	output {
		Array[File] reports = flatten(fastqc.reports)
        File count_matrix_h5ad = SalmonQuantificationCall.count_matrix_h5ad
        File genome_build_json = SalmonQuantificationCall.genome_build_json
        File? raw_count_matrix = SalmonQuantificationCall.raw_count_matrix
	}

}
