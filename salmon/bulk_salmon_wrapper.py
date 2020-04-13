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
    'quant',
    '--index',
    '/opt/grch38_index',
    '--libType',
    'A',
    '-p',
    '{threads}',
    '--output',
    'out'
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

def get_sample_id(filename:str):
    return filename.split("/")[-1].replace("_R1.FASTQ", "")

def rename_file(old_file_name: str, new_file_name: str):
    command = ["mv"]
    command.append(old_file_name)
    command.append(new_file_name)
    check_call(command)

def find_r1_fastq_files(directory: Path) -> Iterable[Path]:
    pattern = '**/*_R1*.{extension}'
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
        r2_fastq_filename = r1_fastq_file.name.replace('_R1', '_R2')
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

    for r1_fastq_file, r2_fastq_file in find_fastq_files(directory):

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
        #tar and zip auxilliary files

        sample_id = get_sample_id(str(r1_fastq_file))

        #Tag output files with sample_id
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
