version development

task FastQC {
    input {
        Directory fastq_dir
        Int threads
    }

    output {
        Directory fastqc_dir = "fastqc_output"
    }

    runtime {
        container: "hubmap/scrna-analysis:latest"
    }

    command {
        /opt/fastqc_wrapper.py ~{fastq_dir} ~{threads}
    }
}