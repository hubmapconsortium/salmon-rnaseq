#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from os import fspath
from subprocess import check_call
from typing import Iterable, Tuple

FASTQ_EXTENSIONS = [
    'fastq',
    'fastq.gz',
]

SALMON_COMMAND = [
    'salmon',
    'alevin',
    '--index',
    '/opt/grch38_index',
    '--libType',
    'A',
    '--output',
    'out',
    '--chromiumV3',
    '--tgMap',
    '/opt/Homo_sapiens.GRCh38.cdna.all.fa.gz.map',
    '-p',
    '{threads}',
]

FOUND_PAIR_COLOR = '\033[01;32m'
UNPAIRED_COLOR = '\033[01;31m'
NO_COLOR = '\033[00m'

def find_r1_fastq_files(directory: Path) -> Iterable[Path]:
    pattern = '**/*_R1_*.{extension}'
    for extension in FASTQ_EXTENSIONS:
        yield from directory.glob(pattern.format(extension=extension))

def find_fastq_files(directory: Path) -> Iterable[Tuple[Path, Path]]:
    """
    Specific to 10X FASTQ filename conventions. Returns all paired R1/R2
    FASTQ files in any subdirectory of 'directory'.

    :param directory:
    :return: Iterable of 2-tuples:
     [0] R1 FASTQ file
     [1] R2 FASTQ file
    """
    for r1_fastq_file in find_r1_fastq_files(directory):
        r2_fastq_filename = r1_fastq_file.name.replace('_R1_', '_R2_')
        r2_fastq_file = r1_fastq_file.with_name(r2_fastq_filename)
        if r2_fastq_file.is_file():
            print(FOUND_PAIR_COLOR + 'Found pair of FASTQ files:' + NO_COLOR)
            print('\t', r1_fastq_file, sep='')
            print('\t', r2_fastq_file, sep='')
            yield r1_fastq_file, r2_fastq_file
        else:
            print(UNPAIRED_COLOR + 'Found unpaired FASTQ file:' + NO_COLOR)
            print('\t', r1_fastq_file, sep='')

def main(threads: int, directory: Path):
    command = [
        piece.format(threads=threads)
        for piece in SALMON_COMMAND
    ]
    for r1_fastq_file, r2_fastq_file in find_fastq_files(directory):
        fastq_extension = [
            '-1',
            fspath(r1_fastq_file),
            '-2',
            fspath(r2_fastq_file),
        ]
        command.extend(fastq_extension)
    print('Running:', ' '.join(command))
    check_call(command)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('-p', '--threads', type=int)
    p.add_argument('directory', type=Path)
    args = p.parse_args()

    main(args.threads, args.directory)
