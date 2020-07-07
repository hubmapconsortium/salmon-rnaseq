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

def rename_prepend_sample_id(path: Path, sample_id: str) -> Path:
    """
    Prepends `sample_id` and '-' to the file/directory name in `path`,
    with and renames the file or directory on disk

    :param path: file/directory to rename
    :param sample_id: new prefix of name of `path`
    :return: the new path
    """
    new_name = f'{sample_id}-{path.name}'
    new_path = path.with_name(new_name)
    return path.rename(new_path)

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

        out_dir = Path('out')
        # inside `out_dir`
        filenames_to_rename = [
            'quant.sf',
            'cmd_info.json',
            'aux_files.tar.gz',
        ]
        for filename in filenames_to_rename:
            rename_prepend_sample_id(out_dir / filename, sample_id)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('-p', '--threads', type=int)
    p.add_argument('directory', type=Path)
    args = p.parse_args()
    print(args.threads, args.directory)

    main(args.threads, args.directory)
