from pathlib import Path
import re
from typing import Iterable, Tuple

FASTQ_R1_PATTERN = re.compile(r'(.*)_(R?)(1)(\.(fq|fastq)(\.gz)?)')

FOUND_PAIR_COLOR = '\033[01;32m'
UNPAIRED_COLOR = '\033[01;31m'
NO_COLOR = '\033[00m'

def get_sample_id_from_r1(file_path: Path) -> str:
    """
    Only supports R1 FASTQ files.

    @param file_path:
    @return:

    >>> get_sample_id_from_r1(Path('whatever/path/to/B001A001_1.fastq'))
    'B001A001'
    >>> get_sample_id_from_r1(Path('whatever/path/to/B001A001_1.fastq.gz'))
    'B001A001'
    >>> get_sample_id_from_r1(Path('whatever/path/to/B001A001_1.fq'))
    'B001A001'
    >>> get_sample_id_from_r1(Path('whatever/path/to/B001A001_1.fq.gz'))
    'B001A001'
    >>> get_sample_id_from_r1(Path('whatever/path/to/B001A001_R1.fastq'))
    'B001A001'
    >>> get_sample_id_from_r1(Path('whatever/path/to/B001A001_R1.fastq.gz'))
    'B001A001'
    >>> get_sample_id_from_r1(Path('whatever/path/to/B001A001_R1.fq'))
    'B001A001'
    >>> get_sample_id_from_r1(Path('whatever/path/to/B001A001_R1.fq.gz'))
    'B001A001'
    >>> get_sample_id_from_r1(Path('whatever/path/to/B001A001_2.fq.gz'))
    Traceback (most recent call last):
      ...
    ValueError: Path did not match R1 FASTQ pattern: whatever/path/to/B001A001_2.fq.gz
    """
    if not FASTQ_R1_PATTERN.match(file_path.name):
        raise ValueError(f'Path did not match R1 FASTQ pattern: {file_path}')
    return FASTQ_R1_PATTERN.sub(r'\1', file_path.name)

def get_r2_fastq(file_path: Path) -> Path:
    """
    @param file_path:
    @return:
    >>> get_r2_fastq(Path('path/to/B001A001_1.fastq'))
    PosixPath('path/to/B001A001_2.fastq')
    >>> get_r2_fastq(Path('path/to/B001A001_1.fastq.gz'))
    PosixPath('path/to/B001A001_2.fastq.gz')
    >>> get_r2_fastq(Path('path/to/B001A001_1.fq'))
    PosixPath('path/to/B001A001_2.fq')
    >>> get_r2_fastq(Path('path/to/B001A001_1.fq.gz'))
    PosixPath('path/to/B001A001_2.fq.gz')
    >>> get_r2_fastq(Path('path/to/B001A001_R1.fastq'))
    PosixPath('path/to/B001A001_R2.fastq')
    >>> get_r2_fastq(Path('path/to/B001A001_R1.fastq.gz'))
    PosixPath('path/to/B001A001_R2.fastq.gz')
    >>> get_r2_fastq(Path('path/to/B001A001_R1.fq'))
    PosixPath('path/to/B001A001_R2.fq')
    >>> get_r2_fastq(Path('path/to/B001A001_R1.fq.gz'))
    PosixPath('path/to/B001A001_R2.fq.gz')
    >>> get_r2_fastq(Path('path/to/B001A001_2.fq.gz'))
    Traceback (most recent call last):
      ...
    ValueError: Path did not match R1 FASTQ pattern: path/to/B001A001_2.fq.gz
    """
    if not FASTQ_R1_PATTERN.match(file_path.name):
        raise ValueError(f'Path did not match R1 FASTQ pattern: {file_path}')
    new_filename = FASTQ_R1_PATTERN.sub(r'\1_\g<2>2\4', file_path.name)
    return file_path.with_name(new_filename)

def find_r1_fastq_files(directory: Path) -> Iterable[Path]:
    for path in directory.glob('**/*'):
        if path.is_file() and FASTQ_R1_PATTERN.match(path.name):
            yield path

def find_fastq_files(directory: Path, verbose=True) -> Iterable[Tuple[Path, Path]]:
    """
    Returns all paired R1/R2 FASTQ files in any subdirectory of 'directory'.

    :param directory:
    :return: Iterable of 2-tuples:
     [0] R1 FASTQ file
     [1] R2 FASTQ file
    """
    for r1_fastq_file in find_r1_fastq_files(directory):
        r2_fastq_file = get_r2_fastq(r1_fastq_file)
        if r2_fastq_file.is_file():
            if verbose:
                print(FOUND_PAIR_COLOR + 'Found pair of FASTQ files:' + NO_COLOR)
                print(f'\t{r1_fastq_file}')
                print(f'\t{r2_fastq_file}')
            yield r1_fastq_file, r2_fastq_file
        else:
            if verbose:
                print(UNPAIRED_COLOR + 'Found unpaired FASTQ file:' + NO_COLOR)
                print(f'\t{r1_fastq_file}')
