#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from os import environ
from subprocess import check_call

SALMON_COMMAND = [
    'salmon',
    'alevin',
    '--index',
    '/opt/gencode.v35.intron-exon.sidx',
    '--libType',
    'A',
    '--output',
    'out',
    '--snareseq',
    '--keepCBFraction',
    '1',
    '--tgMap',
    '/opt/gencode.v35.annotation.expanded.tx2gene.tsv',
    '-p',
    '{threads}',
    '-1',
    '{barcode_umi_fastq}',
    '-2',
    '{transcript_fastq}',
]

def main(threads: int, barcode_umi_fastq: Path, transcript_fastq: Path):
    command = [
        piece.format(
            threads=threads,
            barcode_umi_fastq=barcode_umi_fastq,
            transcript_fastq=transcript_fastq,
        )
        for piece in SALMON_COMMAND
    ]
    print('Running:', ' '.join(command))
    # for Singularity container runtime
    env = environ.copy()
    env['LD_LIBRARY_PATH'] = '/usr/local/lib'
    # /for Singularity container runtime
    check_call(command)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('-p', '--threads', type=int)
    p.add_argument('barcode_umi_fastq', type=Path)
    p.add_argument('transcript_fastq', type=Path)
    args = p.parse_args()

    main(args.threads, args.barcode_umi_fastq, args.transcript_fastq)
