from enum import Enum
from pathlib import Path
from typing import Tuple

ADJ_OUTPUT_DIR = Path('adj_fastq')

# R1 for Salmon
BARCODE_UMI_FASTQ_PATH = Path('barcode_umi.fastq')
# R2 for Salmon
TRANSCRIPT_FASTQ_FILENAME_BASE = 'transcript.fastq'
TRANSCRIPT_FASTQ_PATH = Path(TRANSCRIPT_FASTQ_FILENAME_BASE)
TRANSCRIPT_FASTQ_GZ_PATH = Path(TRANSCRIPT_FASTQ_FILENAME_BASE + '.gz')

def get_adjusted_fastq_paths(directory: Path) -> Tuple[Path, Path]:
    return directory / BARCODE_UMI_FASTQ_PATH, directory / TRANSCRIPT_FASTQ_PATH

class Assay(Enum):
    def __new__(
            cls,
            key: str,
            salmon_option: str,
            barcode_adj_performed: bool,
            keep_all_barcodes: bool,
    ):
        obj = object.__new__(cls)
        obj._value_ = key
        obj.salmon_option = salmon_option
        obj.barcode_adj_performed = barcode_adj_performed
        obj.keep_all_barcodes = keep_all_barcodes
        return obj

    def __str__(self):
        return self.value

    CHROMIUM_V3 = '10x', '--chromiumV3', True, False
    SNARESEQ = 'snareseq', '--snareseq', False, True
    SCISEQ = 'sciseq', '--sciseq', False, True
