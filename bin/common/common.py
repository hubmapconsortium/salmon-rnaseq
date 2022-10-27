from enum import Enum
from pathlib import Path
from subprocess import run
from typing import Tuple

ADJ_OUTPUT_DIR = Path("adj_fastq")

# R1 for Salmon
BARCODE_UMI_FASTQ_PATH = Path("barcode_umi.fastq")
# R2 for Salmon
TRANSCRIPT_FASTQ_FILENAME_BASE = "transcript.fastq"
TRANSCRIPT_FASTQ_PATH = Path(TRANSCRIPT_FASTQ_FILENAME_BASE)
TRANSCRIPT_FASTQ_GZ_PATH = Path(TRANSCRIPT_FASTQ_FILENAME_BASE + ".gz")

# Probably faster than piping through the Python interpreter, even
# though we're reading everything anyway, to write the barcode/UMI FASTQ
GUNZIP_COMMAND = ["gunzip", "-c", "{path}"]


def decompress_fastq(input_path: Path, output_path: Path):
    print("Decompressing", input_path, "to", output_path)
    with open(output_path, "ab") as o:
        command = [piece.format(path=input_path) for piece in GUNZIP_COMMAND]
        run(command, stdout=o, check=True)
    return output_path


def get_adjusted_fastq_paths(directory: Path) -> Tuple[Path, Path]:
    return directory / BARCODE_UMI_FASTQ_PATH, directory / TRANSCRIPT_FASTQ_PATH


class AnnDataLayer(str, Enum):
    SPLICED = "spliced"
    UNSPLICED = "unspliced"
    SPLICED_UNSPLICED_SUM = "spliced_unspliced_sum"


class Assay(Enum):
    def __new__(
        cls,
        key: str,
        salmon_option: str,
        secondary_analysis_layer: AnnDataLayer,
        barcode_adj_performed: bool,
        barcode_adj_r1_r2: bool,
        keep_all_barcodes: bool,
    ):
        obj = object.__new__(cls)
        obj._value_ = key
        obj.salmon_option = salmon_option
        obj.secondary_analysis_layer = secondary_analysis_layer
        obj.barcode_adj_performed = barcode_adj_performed
        obj.barcode_adj_r1_r2 = barcode_adj_r1_r2
        obj.keep_all_barcodes = keep_all_barcodes
        return obj

    def __str__(self):
        return self.value

    CHROMIUM_V2 = (
        "10x_v2",
        "--chromium",
        AnnDataLayer.SPLICED,
        False,
        False,
        False,
    )
    CHROMIUM_V3 = (
        "10x_v3",
        "--chromiumV3",
        AnnDataLayer.SPLICED,
        False,
        False,
        False,
    )
    CHROMIUM_V2_SN = (
        "10x_v2_sn",
        "--chromium",
        AnnDataLayer.SPLICED_UNSPLICED_SUM,
        False,
        False,
        False,
    )
    CHROMIUM_V3_SN = (
        "10x_v3_sn",
        "--chromiumV3",
        AnnDataLayer.SPLICED_UNSPLICED_SUM,
        False,
        False,
        False,
    )
    SNARESEQ = (
        "snareseq",
        "--snareseq",
        AnnDataLayer.SPLICED_UNSPLICED_SUM,
        True,
        False,
        True,
    )
    SCISEQ = (
        "sciseq",
        "--sciseq",
        AnnDataLayer.SPLICED_UNSPLICED_SUM,
        True,
        True,
        True,
    )
    SLIDESEQ = (
        "slideseq",
        "--slideseq",
        AnnDataLayer.SPLICED_UNSPLICED_SUM,
        True,
        False,
        False,
    )
    VISIUM_FFPE = (
        "visium-ffpe",
        "--visium-ffpe",
        AnnDataLayer.SPLICED_UNSPLICED_SUM,
        True,
        False,
        False,
    )