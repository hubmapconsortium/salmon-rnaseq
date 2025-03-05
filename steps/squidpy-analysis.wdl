version 1.0

task SquidPyAnalysis {
    input {
        String assay
        File h5ad_file
        Directory? img_dir
    }

    output {
        File? squidpy_annotated_h5ad = "squidpy_annotated.h5ad"
        File? neighborhood_enrichment_plot = "neighborhood_enrichment.pdf"
        File? co_occurrence_plot = "co_occurrence.pdf"
        File? spatial_plot = "spatial_scatter.pdf"
        File? interaction_matrix_plot = "interaction_matrix.pdf"
        File? centrality_scores_plot = "centrality_scores.pdf"
        File? ripley_plot = "ripley.pdf"
    }

    runtime {
        container: "hubmap/squidpy-analysis:latest"
    }

    command <<<
        /opt/squidpy_entry_point.py ~{assay} ~{h5ad_file} ~{img_dir}
    >>>
}