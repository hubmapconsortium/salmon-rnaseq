#!/usr/bin/env python3
from argparse import ArgumentParser
from itertools import chain
from pathlib import Path
from typing import Iterable

from fastq_utils import Read, fastq_reader, find_grouped_fastq_files

from common import BARCODE_UMI_FASTQ_PATH, TRANSCRIPT_FASTQ_PATH, decompress_fastq

def get_barcode_umi(seq: str) -> str:
    barcode = seq[0:8] + seq[26:32]
    umi = seq[32:41]
    return barcode + umi

def main(fastq_dirs: Iterable[Path], output_dir: Path = Path()):
    buf = output_dir / BARCODE_UMI_FASTQ_PATH
    trf = output_dir / TRANSCRIPT_FASTQ_PATH

    all_fastqs = chain.from_iterable(
        find_grouped_fastq_files(fastq_dir, 2)
        for fastq_dir in fastq_dirs
    )

    with open(buf, 'w') as cbo:
        for transcript_fastq, barcode_umi_fastq in all_fastqs:
            decompress_fastq(transcript_fastq, trf)
            print('Extracting barcodes from', barcode_umi_fastq)
            barcode_umi_reader = fastq_reader(barcode_umi_fastq)
            for br in barcode_umi_reader:
                new_br = Read(
                    read_id=br.read_id,
                    seq=get_barcode_umi(br.seq),
                    unused=br.unused,
                    qual=get_barcode_umi(br.qual),
                )
                print(new_br.serialize(), file=cbo)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('directory', type=Path)
    args = p.parse_args()

    main(args.directory)
