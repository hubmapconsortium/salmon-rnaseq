#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from os import environ, fspath
from subprocess import check_call

from fastq_utils import find_fastq_files

FASTQ_EXTENSIONS = [
    'fastq',
    'fastq.gz',
    'fq.gz',
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
    env = environ.copy()
    env['LD_LIBRARY_PATH'] = '/usr/local/lib'
    check_call(command)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('-p', '--threads', type=int)
    p.add_argument('directory', type=Path)
    args = p.parse_args()

    main(args.threads, args.directory)
