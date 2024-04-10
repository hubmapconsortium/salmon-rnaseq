## Use double '#' for workflow-level comments
## This workflow implements a one-task workflow

# write the WDL version number 'version 1.0' -- 1
# possible to write 'WDL developent' as a version number as well
version development

# create a workflow named 'HelloWorld' -- 2
import "./steps/salmon-quantification.wdl" as SalmonQuantification
import "./steps/fastqc.wdl" as FastQC
import "./steps/scanpy-analysis.wdl" as ScanPyAnalysis
import "./steps/scvelo-analysis.wdl" as ScVeloAnalysis
import "./steps/squidpy-analysis.wdl" as SquidPyAnalysis
import "./steps/compute-qc-metrics.wdl" as ComputeQCMetrics

workflow SalmonRNAseq {
    input {
        Array[Directory] fastq_dir
        Directory? img_dir
        Directory? metadata_dir
        String assay
        Int threads
        Int? expected_cell_count
        Boolean? keep_all_barcodes
    }

    call SalmonQuantification.SalmonQuantification as SalmonQuantificationCall {
        input:
            fastq_dir = fastq_dir,
            img_dir = img_dir,
            metadata_dir = metadata_dir,
            assay = assay,
            threads = threads,
            expected_cell_count = expected_cell_count,
            keep_all_barcodes = keep_all_barcodes
    }

    scatter (fastq in fastq_dir) {
        call FastQC.FastQC as FastQCCall {
            input:
                fastq_dir = fastq,
                threads = threads
        }
    }

    call ScanPyAnalysis.ScanPyAnalysis as ScanPyAnalysisCall {
        input:
            assay = assay,
            h5ad_file = SalmonQuantificationCall.count_matrix_h5ad
    }

    call ScVeloAnalysis.ScVeloAnalysis as ScVeloAnalysisCall {
        input:
            spliced_h5ad_file = SalmonQuantificationCall.count_matrix_h5ad,
            assay_name = assay
    }

    call SquidPyAnalysis.SquidPyAnalysis as SquidPyAnalysisCall {
        input:
            assay = assay,
            h5ad_file = ScanPyAnalysisCall.filtered_data_h5ad,
            img_dir = img_dir
    }

    call ComputeQCMetrics.ComputeQCMetrics as ComputeQCMetricsCall {
        input:
            assay = assay,
            h5ad_primary = SalmonQuantificationCall.count_matrix_h5ad,
            h5ad_secondary = ScanPyAnalysisCall.filtered_data_h5ad,
            salmon_dir = SalmonQuantificationCall.salmon_output
    }
}