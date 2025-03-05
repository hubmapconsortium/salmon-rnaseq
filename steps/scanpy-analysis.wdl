version 1.0

task ScanPyAnalysis {
    input {
        String assay
        File h5ad_file
    }

    output {
        File filtered_data_h5ad = "secondary_analysis.h5ad"
        File dispersion_plot = "dispersion_plot.pdf"
        File umap_plot = "umap_by_leiden_cluster.pdf"
        File? spatial_plot = "spatial_pos_by_leiden_cluster.pdf"
        File umap_density_plot = "umap_embedding_density.pdf"
        File marker_gene_plot_t_test = "marker_genes_by_cluster_t_test.pdf"
        File marker_gene_plot_logreg = "marker_genes_by_cluster_logreg.pdf"
    }

    runtime {
        container: "hubmap/scrna-analysis:latest"
    }

    command <<<
        /opt/scanpy_entry_point.py ~{assay} ~{h5ad_file}
    >>>
}