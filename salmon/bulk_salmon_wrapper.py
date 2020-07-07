#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from subprocess import check_call

from fastq_utils import (
    get_sample_id_from_r1,
    find_fastq_files,
)

SALMON_COMMAND = [
    'salmon',
    'quant',
    '--index',
    '/opt/grch38_index',
    '--libType',
    'A',
    '-p',
    '{threads}',
    '--output',
    'out',
    '-1',
    '{fastq_r1}',
    '-2',
    '{fastq_r2}',
]

TAR_AND_ZIP_COMMAND = [
    'tar',
    '-czvf',
    'out/aux_files.tar.gz',
    'out/aux_info',
]

FOUND_PAIR_COLOR = '\033[01;32m'
UNPAIRED_COLOR = '\033[01;31m'
NO_COLOR = '\033[00m'

def rename_file(old_file_name: str, new_file_name: str):
    command = ["mv"]
    command.append(old_file_name)
    command.append(new_file_name)
    check_call(command)

def main(threads: int, directory: Path):
    for r1_fastq_file, r2_fastq_file in find_fastq_files(directory):
        command = [
            piece.format(
                threads=threads,
                fastq_r1=r1_fastq_file,
                fastq_r2=r2_fastq_file,
            )
            for piece in SALMON_COMMAND
        ]
        print('Running:', command)
        check_call(command)

        check_call(TAR_AND_ZIP_COMMAND)
        # tar and zip auxilliary files

        sample_id = get_sample_id_from_r1(r1_fastq_file)

        # Tag output files with sample_id
        rename_file("out/quant.sf", "out/" + sample_id + "-quant.sf")
        rename_file("out/cmd_info.json", "out/" + sample_id + "-cmd_info.json")
        rename_file("out/aux_files.tar.gz", "out/" + sample_id + "-aux_files.tar.gz")

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('-p', '--threads', type=int)
    p.add_argument('directory', type=Path)
    args = p.parse_args()
    print(args.threads, args.directory)

    main(args.threads, args.directory)
