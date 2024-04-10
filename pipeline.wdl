## Use double '#' for workflow-level comments
## This workflow implements a one-task workflow

# write the WDL version number 'version 1.0' -- 1
# possible to write 'WDL developent' as a version number as well
version development

# create a workflow named 'HelloWorld' -- 2
import "./steps/salmon-quantification.wdl" as SalmonQuantification
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

    call SalmonQuantification.SalmonQuantification {
        input:
            fastq_dir = fastq_dir,
            img_dir = img_dir,
            metadata_dir = metadata_dir,
            assay = assay,
            threads = threads,
            expected_cell_count = expected_cell_count,
            keep_all_barcodes = keep_all_barcodes
    }
}
