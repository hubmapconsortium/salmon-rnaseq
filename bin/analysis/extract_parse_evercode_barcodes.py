#!/usr/bin/env python3
from argparse import ArgumentParser
from itertools import chain
from pathlib import Path
from typing import Iterable

import manhole
from fastq_utils import Read, fastq_reader, find_grouped_fastq_files

from common import BARCODE_UMI_FASTQ_PATH, TRANSCRIPT_FASTQ_PATH

BARCODE_LENGTH = 8
# order as per vendor documentation
BARCODE_STARTS = [78, 48, 10]
BARCODE_SEGMENTS = [slice(start, start + BARCODE_LENGTH) for start in BARCODE_STARTS]
UMI_SEGMENT = slice(0, 10)

BARCODE_QUAL_DUMMY = "F" * BARCODE_LENGTH * len(BARCODE_STARTS)


def main(
    fastq_dirs: Iterable[Path],
    output_dir: Path = Path(),
):
    buf = output_dir / BARCODE_UMI_FASTQ_PATH
    trf = output_dir / TRANSCRIPT_FASTQ_PATH

    all_fastqs = chain.from_iterable(
        find_grouped_fastq_files(fastq_dir, 2) for fastq_dir in fastq_dirs
    )

    with open(buf, "w") as cbo, open(trf, "w") as tro:
        for transcript_fastq, barcode_umi_fastq in all_fastqs:
            usable_count = 0
            i = 0
            print("Correcting barcodes in", transcript_fastq, "and", barcode_umi_fastq)
            transcript_reader = fastq_reader(transcript_fastq)
            barcode_umi_reader = fastq_reader(barcode_umi_fastq)
            for i, (tr, br) in enumerate(zip(transcript_reader, barcode_umi_reader), 1):
                barcode_pieces = [br.seq[s] for s in BARCODE_SEGMENTS]
                bc_qual_pieces = [br.qual[s] for s in BARCODE_SEGMENTS]
                usable_count += 1
                umi_seq = br.seq[UMI_SEGMENT]
                umi_qual = br.qual[UMI_SEGMENT]
                new_seq = "".join(barcode_pieces + [umi_seq])
                new_qual = "".join(bc_qual_pieces + [umi_qual])
                new_br = Read(
                    read_id=br.read_id,
                    seq=new_seq,
                    unused=br.unused,
                    qual=new_qual,
                )
                print(tr.serialize(), file=tro)
                print(new_br.serialize(), file=cbo)

            print("Total count:", i)
            print("Usable count:", usable_count)
            print("Proportion:", usable_count / i)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("fastq_dirs", type=Path, nargs="+")
    args = p.parse_args()

    main(fastq_dirs=args.fastq_dirs)
