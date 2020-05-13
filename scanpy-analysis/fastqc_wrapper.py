#!/usr/bin/env python3
from argparse import ArgumentParser
from collections import defaultdict
from os import fspath, walk
from pathlib import Path
from subprocess import check_call
from typing import Dict, List, Iterable

FASTQ_PATTERNS = [
    '*.fastq',
    '*.fastq.gz',
]
FASTQC_COMMAND_TEMPLATE = [
    'fastqc',
    '--outdir',
    '{out_dir}',
]

def find_fastq_files(directory: Path) -> Iterable[Path]:
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            if any(filepath.match(pattern) for pattern in FASTQ_PATTERNS):
                yield filepath.relative_to(directory)

def collect_fastq_files_by_directory(directory: Path) -> Dict[Path, List[Path]]:
    """
    Walk `directory`, finding all FASTQ files. Group these by the directory
    the files are in, so we can create the same directory structure for the
    output of FastQC.

    :param directory: Path to directory containing FASTQ files
    :return: Mapping of *relative* directory names to lists of absolute
      FASTQ file paths. The relative directory names are used to stage
      output directories matching the same structure as the input directory
      tree; the FASTQ paths are passed directly to FastQC
    """
    files_by_directory = defaultdict(list)
    for fastq_file in find_fastq_files(directory):
        # fastq_file is relative to `directory`, make absolute again
        files_by_directory[fastq_file.parent].append(directory / fastq_file)
    return files_by_directory

def main(directory: Path):
    fastq_files_by_directory = collect_fastq_files_by_directory(directory)
    print('Found', len(fastq_files_by_directory), 'directories containing FASTQ files')

    fastqc_out_dir = Path('fastqc_output')
    for directory, files in fastq_files_by_directory.items():
        subdir = fastqc_out_dir / directory
        subdir.mkdir(exist_ok=True, parents=True)
        command = [
            piece.format(out_dir=subdir)
            for piece in FASTQC_COMMAND_TEMPLATE
        ]
        command.extend(fspath(fastq_file) for fastq_file in files)
        print('Running', ' '.join(command))
        check_call(command)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('directory', type=Path)
    args = p.parse_args()

    main(args.directory)
