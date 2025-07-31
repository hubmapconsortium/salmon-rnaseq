version development

task ScVeloAnalysis {
    input {
        File spliced_h5ad_file
        String assay_name
    }

    output {
        File? annotated_h5ad_file = "scvelo_annotated.h5ad"
        File? embedding_grid_plot = "scvelo_embedding_grid.pdf"
    }

    runtime {
        container: "hubmap/scrna-analysis:latest"
    }

    command {
        /opt/scvelo_analysis.py ~{spliced_h5ad_file} ~{assay_name}
    }
}