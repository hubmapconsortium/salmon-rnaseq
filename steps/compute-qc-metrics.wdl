version development

task ComputeQCMetrics {
    input {
        String assay
        File h5ad_primary
        File h5ad_secondary
        Directory salmon_dir
    }

    output {
        File scanpy_qc_results = "qc_results.hdf5"
        File qc_metrics = "qc_results.json"
    }

    runtime {
        container: "hubmap/scrna-analysis:latest"
    }

    command {
        /opt/compute_qc_metrics.py ~{assay} ~{h5ad_primary} ~{h5ad_secondary} ~{salmon_dir}
    }
}