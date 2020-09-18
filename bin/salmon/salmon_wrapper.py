#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from os import environ
from subprocess import check_call
from typing import Iterable, List, Sequence, Tuple

from fastq_utils import find_fastq_files

from common import (
    Assay,
    BARCODE_UMI_FASTQ_PATH,
    TRANSCRIPT_FASTQ_GZ_PATH,
    TRANSCRIPT_FASTQ_PATH,
)

SALMON_COMMAND = [
    'salmon',
    'alevin',
    '--index',
    '/opt/gencode.v35.intron-exon.sidx',
    '--libType',
    'A',
    '--output',
    'salmon_out',
    '--dumpMtx',
    '{salmon_option}',
    '--tgMap',
    '/opt/gencode.v35.annotation.expanded.tx2gene.tsv',
    '-p',
    '{threads}',
]

def find_adj_fastq_files(directory: Path) -> Tuple[Path, Path]:
    # not general enough to implement in fastq-utils
    barcode_umi_fastq = directory / BARCODE_UMI_FASTQ_PATH

    transcript_fastq = directory / TRANSCRIPT_FASTQ_PATH
    transcript_fastq_gz = directory / TRANSCRIPT_FASTQ_GZ_PATH

    if transcript_fastq.is_file():
        return barcode_umi_fastq, transcript_fastq
    elif transcript_fastq_gz.is_file():
        return barcode_umi_fastq, transcript_fastq_gz
    else:
        message = f"Couldn't find {TRANSCRIPT_FASTQ_PATH} or {TRANSCRIPT_FASTQ_GZ_PATH} in {directory}"
        raise ValueError(message)

def main(assay: Assay, orig_fastq_dir: Iterable[Path], adj_fastq_dir: Path, threads: int):
    command = [
        piece.format(
            salmon_option=assay.salmon_option,
            threads=threads,
        )
        for piece in SALMON_COMMAND
    ]

    fastq_pairs: List[Sequence[Path]]
    if assay.barcode_adj_performed:
        if assay.barcode_adj_r1_r2:
            fastq_pairs = list(find_fastq_files([adj_fastq_dir], 2))
        else:
            fastq_pairs = [find_adj_fastq_files(adj_fastq_dir)]
    else:
        fastq_pairs = list(find_fastq_files(orig_fastq_dir, 2))

    if assay.keep_all_barcodes:
        command.extend(['--keepCBFraction', '1'])

    for r1_fastq_file, r2_fastq_file in fastq_pairs:
        fastq_extension = [
            '-1',
            r1_fastq_file,
            '-2',
            r2_fastq_file,
        ]
        command.extend(fastq_extension)

    print('Running:', ' '.join(str(x) for x in command))
    env = environ.copy()
    # Necessary for Singularity; this environment variable isn't
    # set by that container runtime but is required to run Salmon
    env['LD_LIBRARY_PATH'] = '/usr/local/lib'
    check_call(command)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('assay', choices=list(Assay), type=Assay)
    p.add_argument('adj_fastq_dir', type=Path)
    p.add_argument('orig_fastq_dir', type=Path, nargs='+')
    p.add_argument('-p', '--threads', type=int)
    args = p.parse_args()

    main(args.assay, args.orig_fastq_dir, args.adj_fastq_dir, args.threads)
