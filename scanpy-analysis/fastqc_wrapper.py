#!/usr/bin/env python3

from pathlib import Path
from os import fspath
from argparse import ArgumentParser
from subprocess import check_call
from typing import Iterable, Tuple

FASTQ_EXTENSIONS = [
    'fastq',
    'fastq.gz',
]

def find_fastq_files(directory: Path) -> Iterable[Path]:
    pattern = '*.{extension}'
    for extension in FASTQ_EXTENSIONS:
        yield from directory.glob(pattern.format(extension=extension))

def main(directory: Path):
    command = ["fastqc", "--outdir", "."]
    for fastq_file in find_fastq_files(directory):
        command.append(str(fastq_file))

    print(command)

    print('Running:', ' '.join(command))
    check_call(command)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('directory', type=Path)
    args = p.parse_args()

    main(args.directory)
