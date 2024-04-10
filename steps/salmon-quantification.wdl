version development

workflow SalmonQuantification {
    input {
        Array[Directory] fastq_dir
        Directory? img_dir
        Directory? metadata_dir
        String assay
        Int threads
        Int? expected_cell_count
        Boolean? keep_all_barcodes
    }

    output {
        Directory salmon_output = Salmon.output_dir
        File count_matrix_h5ad = AnnotateCells.annotated_h5ad_file
        File? raw_count_matrix = AlevinToAnndata.raw_expr_h5ad
        File genome_build_json = AlevinToAnndata.genome_build_json
    }

    call AdjustBarcodes{
      input:
        assay = assay,
        fastq_dir = fastq_dir
    }

    call TrimReads {
        input:
            assay = assay,
            adj_fastq_dir = AdjustBarcodes.adj_fastq_dir,
            orig_fastq_dirs = fastq_dir,
            threads = threads
    }

    call Salmon {
        input:
            orig_fastq_dirs = fastq_dir,
            trimmed_fastq_dir = TrimReads.trimmed_fastq_dir,
            assay = assay,
            threads = threads,
            expected_cell_count = expected_cell_count,
            keep_all_barcodes = keep_all_barcodes
    }

    call AlevinToAnndata {
        input:
            assay = assay,
            alevin_dir = Salmon.output_dir
    }

    call AnnotateCells {
        input:
            assay = assay,
            orig_fastq_dirs = fastq_dir,
            h5ad_file = AlevinToAnndata.expr_h5ad,
            img_dir = img_dir,
            metadata_dir = metadata_dir,
            metadata_json = AdjustBarcodes.metadata_json
    }
}

task AdjustBarcodes {
    input {
        String assay
        Array[Directory] fastq_dir
    }

    output {
        Directory adj_fastq_dir = "adj_fastq"
        File? metadata_json = "metadata.json"
    }

    command {
        /opt/adjust_barcodes.py ~{assay} directory ~{sep(" ", fastq_dir)}
    }

    runtime {
        container: "hubmap/scrna-barcode-adj:latest"
    }
}

task TrimReads {
    input {
        String assay
        Directory adj_fastq_dir
        Array[Directory] orig_fastq_dirs
        Int threads
    }

    output {
        Directory trimmed_fastq_dir = "trimmed"
    }

    runtime {
        container: "hubmap/scrna-trim-reads:latest"
    }

    command {
        /opt/trim_reads.py ~{assay} ~{adj_fastq_dir} ~{sep(" ", orig_fastq_dirs)}
    }
}

task Salmon {
    input {
        String assay
        Directory trimmed_fastq_dir
        Array[Directory] orig_fastq_dirs
        Int threads
        Int? expected_cell_count
        Boolean? keep_all_barcodes
    }

    output {
        Directory output_dir = "salmon_out"
    }

    runtime {
        container: "hubmap/salmon-grch38:latest"
    }

    command {
        /opt/salmon_wrapper.py ~{assay} ~{trimmed_fastq_dir} ~{sep(" ", orig_fastq_dirs)} --threads ~{threads} ~{if defined(expected_cell_count) then "--expected-cell-count " + expected_cell_count else ""} ~{if defined(keep_all_barcodes) then "--keep-all-barcodes " + keep_all_barcodes else ""}
    }
}

task AlevinToAnndata {
    input {
        String assay
        Directory alevin_dir
    }

    output {
        File? raw_expr_h5ad = "raw_expr.h5ad"
        File expr_h5ad = "expr.h5ad"
        File genome_build_json = "genome_build.json"
    }

    runtime {
        container: "hubmap/scrna-analysis:latest"
    }

    # Need to fix this command
    command {
        /opt/alevin_to_anndata.py ~{assay} ~{alevin_dir}
    }
}

task AnnotateCells {
    input {
        String assay
        File h5ad_file
        Array[Directory] orig_fastq_dirs
        Directory? img_dir
        Directory? metadata_dir
        File? metadata_json
    }

    output {
        File annotated_h5ad_file = "expr.h5ad"
    }

    runtime {
        container: "hubmap/scrna-analysis:latest"
    }

    # Need to fix this command
    command {
        /opt/annotate_cells.py
    }
}