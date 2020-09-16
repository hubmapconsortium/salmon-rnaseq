#!/usr/bin/env python3
from argparse import ArgumentParser
from os import fspath
from pathlib import Path
from subprocess import check_call

from fastq_utils import find_fastq_files, get_sample_id_from_r1

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
]

TAR_AND_ZIP_COMMAND = [
    'tar',
    '-czvf',
    'out/aux_files.tar.gz',
    'out/aux_info',
]

def rename_file(old_file_name: str, new_file_name: str):
    command = ["mv"]
    command.append(old_file_name)
    command.append(new_file_name)
    check_call(command)

def main(threads: int, directory: Path):
    for r1_fastq_file, r2_fastq_file in find_fastq_files([directory], 2):

        command = [
            piece.format(threads=threads)
            for piece in SALMON_COMMAND
        ]

        fastq_extension = [
            '-1',
            fspath(r1_fastq_file),
            '-2',
            fspath(r2_fastq_file),
        ]

        command.extend(fastq_extension)
        print('Running:', command)
        check_call(command)

        check_call(TAR_AND_ZIP_COMMAND)
        #tar and zip auxiliary files

        sample_id = get_sample_id_from_r1(r1_fastq_file)

        #Tag output files with sample_id
        rename_file("out/quant.sf", f"out/{sample_id}-quant.sf")
        rename_file("out/cmd_info.json", f"out/{sample_id}-cmd_info.json")
        rename_file("out/aux_files.tar.gz", f"out/{sample_id}-aux_files.tar.gz")

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('-p', '--threads', type=int)
    p.add_argument('directory', type=Path)
    args = p.parse_args()
    print(args.threads, args.directory)

    main(args.threads, args.directory)
