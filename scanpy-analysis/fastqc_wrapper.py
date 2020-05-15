#!/usr/bin/env python3
from argparse import ArgumentParser
from collections import defaultdict
from os import fspath, walk
from pathlib import Path
from subprocess import check_call
from typing import Dict, List, Tuple, Iterable
from multiprocessing import Pool

FASTQ_PATTERNS = [
    '*.fastq',
    '*.fastq.gz',
    '*.fq.gz',

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

def single_file_fastqc(fastq_file_and_subdir: Tuple[Path, Path]) -> str:
    #Run fastqc on a single fastq file
    #Takes an absolute path to the input file and a relative path to the output subdirectory
    #Returns a str because imap_unordered seemse to really want this function to be fruitful

    command = [
        piece.format(out_dir=fastq_file_and_subdir[1])
        for piece in FASTQC_COMMAND_TEMPLATE
    ]
    command.append(str(fastq_file_and_subdir[0]))
    print('Running', ' '.join(command))
    check_call(command)
    return command + " completed"

def main(directory: Path, threads:int):
    #Crawl directory, create appropriate output subdirectories based on input directory structure
    #Append output files to a list to pass to Pool.imap_unordered
    fastq_files_by_directory = collect_fastq_files_by_directory(directory)
    print('Found', len(fastq_files_by_directory), 'directories containing FASTQ files')

    fastqc_out_dir = Path('fastqc_output')

    fastq_files_and_subdirs = []

    for directory, files in fastq_files_by_directory.items():
        subdir = fastqc_out_dir / directory
        subdir.mkdir(exist_ok=True, parents=True)
        for file in files:
            fastq_files_and_subdirs.append((file, subdir))

    with Pool(processes=threads) as p:
        for message in p.imap_unordered(single_file_fastqc, fastq_files_and_subdirs):
            print(message)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('directory', type=Path)
    p.add_argument('threads', type=int)
    args = p.parse_args()

    main(args.directory, args.threads)
