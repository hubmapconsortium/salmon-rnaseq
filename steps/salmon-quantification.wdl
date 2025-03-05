version 1.0

workflow salmon_quantification {
    input {
        Int threads
        Int mem_gb
        File fastq1
        File fastq2
        String assay
        String? run_id
        String? protocol
        String? species
        String? img_dir
        String? metadata_dir
        String? organism
        Boolean? keep_all_barcodes
        Int? expected_cell_count
    }

    # Step 1: Adjust Barcodes
    call adjust_barcodes {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            orig_fastqss = sub(fastq1, "/[^/]+$", ""),
            assay = assay
    }

    # Step 2: Trim Reads
    call trim_reads {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            adj_fastqs = adjust_barcodes.adj_fastqs,
            assay = assay
    }

    # Step 3: Salmon Quantification
    call salmon {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            trimmed_fastqs = trim_reads.trimmed_fastqs,
            assay = assay,
            threads = threads,
            mem_gb = mem_gb,
            expected_cell_count = expected_cell_count,
            keep_all_barcodes = keep_all_barcodes,
            organism = organism
    }

    # Step 4: Convert Alevin to AnnData
    call alevin_to_anndata {
        input:
            assay = assay,
            alevin_inputs = salmon.alevin_outputs,
            species = species
            
    }

    # Step 5: Annotate Cells
    call annotate_cells {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            assay = assay,
            h5ad_file = alevin_to_anndata.expr_h5ad,
            img_dir = select_first([img_dir, "default_value"]),
            metadata_dir = select_first([metadata_dir, "default_value"]),
            metadata_json = select_first([adjust_barcodes.metadata_json, "default_value"])
    }

    output {
        File count_matrix_h5ad = alevin_to_anndata.expr_h5ad
        File? raw_count_matrix = alevin_to_anndata.raw_expr_h5ad
        File genome_build_json = alevin_to_anndata.genome_build_json
        File? annotated_h5ad = annotate_cells.annotated_h5ad_file
    }
}

# Task Definitions

task adjust_barcodes {
    input {
        File fastq1
        File fastq2
        String orig_fastqss
        String assay
    }
    
    command {
        # Command for barcode adjustment
        /opt/adjust_barcodes.py ~{assay} ~{orig_fastqss}
    }

    output {
        String adj_fastqs = "/adj_fastq"
        File? metadata_json = "metadata.json"
    }

    runtime {
        docker: "hubmap/scrna-barcode-adj:latest"
    }
}

task trim_reads {
    input {
        File fastq1
        File fastq2
        String adj_fastqs
        String assay
    }

    command {
        # Command for trimming reads
        mkdir -p "adj_fastqs"
        mkdir -p "orig_fastqss"
        mv ~{fastq1} "./orig_fastqss/"
        mv ~{fastq2} "./orig_fastqss/"
        /opt/trim_reads.py ~{assay} "./adj_fastqs" "./orig_fastqss"
    }

    output {
        # Capture all trimmed fastq files
        Array[File] trimmed_fastqs = glob("trimmed/1/*.gz")
        File? metadata_json = "metadata.json"
    }

    runtime {
        docker: "hubmap/scrna-trim-reads:latest"
    }
}

task salmon {
    input {
        Array[File] trimmed_fastqs
        File fastq1
        File fastq2
        String assay
        Int threads
        Int mem_gb
        String? expected_cell_count
        String? keep_all_barcodes
        String? organism
    }

    command {
        mkdir -p "trimmed_fastqs"
        for fastq in ~{sep=' ' trimmed_fastqs}; do
            mv "$fastq" "./trimmed_fastqs/"
        done
        mkdir -p "orig_fastqss"
        mv ~{fastq1} "./orig_fastqss/"
        mv ~{fastq2} "./orig_fastqss/"
        # Command for Salmon quantification (human)
        /opt/salmon_wrapper.py ~{assay} "./trimmed_fastqs/" "./orig_fastqss/" --threads ~{threads} ~{if defined(expected_cell_count) then "--expected-cell-count " + expected_cell_count else ""} ~{if defined(keep_all_barcodes) then "--keep-all-barcodes " + keep_all_barcodes else ""}
    }

    output {
        Array[File] alevin_outputs = glob("salmon_out/alevin/*")
    }

    runtime {
        docker: "hubmap/salmon-grch38:latest"
        memory: "28G"
        cpu: "4"
    }
}

task alevin_to_anndata {
    input {
        String assay
        String? species
        Array[File] alevin_inputs
    }

    command {
        mkdir -p "alevin"
        for file in ~{sep=' ' alevin_inputs}; do
            mv "$file" "./alevin/"
        done
        # Command to convert Alevin output to AnnData object, image default read data from ./alevin
        /opt/alevin_to_anndata.py ~{assay} "./"
    }

    output {
        File? raw_expr_h5ad = "raw_expr.h5ad"
        File expr_h5ad = "expr.h5ad"
        File genome_build_json = "genome_build.json"
    }

    runtime {
        docker: "hubmap/scrna-analysis:latest"
    }
}


task annotate_cells {
    input {
        File fastq1
        File fastq2
        String assay
        File h5ad_file
        String img_dir
        String metadata_dir
        String metadata_json
    }

    command {
        mkdir -p "orig_fastqss"
        mv ~{fastq1} "./orig_fastqss/"
        mv ~{fastq2} "./orig_fastqss/"
        # Command to annotate cells in the AnnData file
        /opt/annotate_cells.py ~{assay} ~{h5ad_file} "./orig_fastqss/" ~{if defined(img_dir) then "--img_dir " + img_dir else ""} ~{if defined(metadata_dir) then "--metadata_dir " + metadata_dir else ""} ~{if defined(metadata_json) then "--metadata_json " + metadata_json else ""}
    }

    output {
        File annotated_h5ad_file = "expr.h5ad"
    }

    runtime {
        docker: "hubmap/scrna-analysis:latest"
    }
}