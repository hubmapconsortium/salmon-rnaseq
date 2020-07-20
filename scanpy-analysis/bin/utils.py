from enum import Enum
import bz2
import gzip
import lzma
from os import PathLike
from pathlib import Path
from typing import Iterable, Tuple

class FileType(Enum):
    def __new__(cls, filetype, open_function):
        obj = object.__new__(cls)
        obj._value_ = filetype
        obj.open_function = open_function
        return obj

    GZ = ("gz", gzip.open)
    BZ2 = ("bz2", bz2.open)
    XZ = ("xz", lzma.open)
    TEXT = ("txt", open)

def get_file_type_by_extension(file_path: Path) -> FileType:
    suffix = file_path.suffix.lstrip(".")
    try:
        return FileType(suffix)
    except ValueError:
        # No special suffix, assume text
        return FileType.TEXT

def smart_open(file_path: PathLike, mode="rt", *args, **kwargs):
    file_type = get_file_type_by_extension(Path(file_path))
    return file_type.open_function(file_path, mode, *args, **kwargs)

FASTQ_EXTENSIONS = [
    'fastq',
    'fastq.gz',
    'fq',
    'fq.gz',
]

FOUND_PAIR_COLOR = '\033[01;32m'
UNPAIRED_COLOR = '\033[01;31m'
NO_COLOR = '\033[00m'

def find_r1_fastq_files(directory: Path) -> Iterable[Path]:
    pattern = '**/*_R1*.{extension}'
    for extension in FASTQ_EXTENSIONS:
        yield from directory.glob(pattern.format(extension=extension))

def find_fastq_files(directory: Path) -> Iterable[Tuple[Path, Path]]:
    """
    :param directory:
    :return: Iterable of 2-tuples:
     [0] R1 FASTQ file
     [1] R2 FASTQ file
    """
    for r1_fastq_file in find_r1_fastq_files(directory):
        r2_fastq_filename = r1_fastq_file.name.replace('_R1', '_R2')
        r2_fastq_file = r1_fastq_file.with_name(r2_fastq_filename)

        if r2_fastq_file.is_file():
            print(FOUND_PAIR_COLOR + 'Found pair of FASTQ files:' + NO_COLOR)
            print(f'\t{r1_fastq_file}')
            print(f'\t{r2_fastq_file}')
            yield r1_fastq_file, r2_fastq_file
        else:
            print(UNPAIRED_COLOR + 'Found unpaired FASTQ file:' + NO_COLOR)
            print(f'\t{r1_fastq_file}')
