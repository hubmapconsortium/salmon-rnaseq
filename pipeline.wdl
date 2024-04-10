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

    output {
        Directory salmon_output = SalmonQuantificationCall.salmon_output
        File count_matrix_h5ad = SalmonQuantificationCall.count_matrix_h5ad
        File? raw_count_matrix = SalmonQuantificationCall.raw_count_matrix
        File genome_build_json = SalmonQuantificationCall.genome_build_json
        Array[Directory] fastqc_dir = FastQCCall.fastqc_dir
        File scanpy_qc_results = ComputeQCMetricsCall.scanpy_qc_results
        File qc_report = ComputeQCMetricsCall.qc_metrics
        File dispersion_plot = ScanPyAnalysisCall.dispersion_plot
        File umap_plot = ScanPyAnalysisCall.umap_plot
        File umap_density_plot = ScanPyAnalysisCall.umap_density_plot
        File? spatial_plot = ScanPyAnalysisCall.spatial_plot
        File filtered_data_h5ad = ScanPyAnalysisCall.filtered_data_h5ad
        File marker_gene_plot_t_test = ScanPyAnalysisCall.marker_gene_plot_t_test
        File marker_gene_plot_logreg = ScanPyAnalysisCall.marker_gene_plot_logreg
        File? scvelo_annotated_h5ad = ScVeloAnalysisCall.annotated_h5ad_file
        File? scvelo_embedding_grid_plot = ScVeloAnalysisCall.embedding_grid_plot
        File? squidpy_annotated_h5ad = SquidPyAnalysisCall.squidpy_annotated_h5ad
        File? neighborhood_enrichment_plot = SquidPyAnalysisCall.neighborhood_enrichment_plot
        File? co_occurrence_plot = SquidPyAnalysisCall.co_occurrence_plot
        File? interaction_matrix_plot = SquidPyAnalysisCall.interaction_matrix_plot
        File? centrality_scores_plot = SquidPyAnalysisCall.centrality_scores_plot
        File? ripley_plot = SquidPyAnalysisCall.ripley_plot
        File? squidpy_spatial_plot = SquidPyAnalysisCall.spatial_plot
    }
}