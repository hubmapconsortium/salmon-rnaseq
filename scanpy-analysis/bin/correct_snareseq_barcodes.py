#!/usr/bin/env python3
from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Set

import barcodeutils as bu

from utils import find_fastq_files, smart_open

BARCODE_LENGTH = 8
BARCODE_STARTS = [10, 48, 86]
BARCODE_SEGMENTS = [slice(start, start + BARCODE_LENGTH) for start in BARCODE_STARTS]
UMI_SEGMENT = slice(0, 10)

BARCODE_QUAL_DUMMY = 'F' * BARCODE_LENGTH * len(BARCODE_STARTS)

@dataclass
class Read:
    read_id: str
    seq: str
    unused: str
    qual: str

    def serialize(self):
        return '\n'.join([self.read_id, self.seq, self.unused, self.qual])

revcomp_table = str.maketrans("ACTG", "TGAC")

def revcomp(seq: str) -> str:
    return seq.translate(revcomp_table)[::-1]

def fastq_reader(fastq_file: Path) -> Iterable[Read]:
    with smart_open(fastq_file) as f:
        while True:
            lines = [f.readline().strip() for _ in range(4)]
            if not all(lines):
                return
            yield Read(*lines)

def read_barcode_allowlist(barcode_filename: Path) -> Set[str]:
    with open(barcode_filename) as f:
        return set(f.read().split())

def main(fastq_dirs: Iterable[Path], barcode_filename: Path):
    barcode_allowlist = read_barcode_allowlist(barcode_filename)
    correcter = bu.BarcodeCorrecter(barcode_allowlist, edit_distance=2)

    with open('barcode_umi.fastq', 'w') as cbo, open('transcript.fastq', 'w') as tro:
        for transcript_fastq, barcode_umi_fastq in find_fastq_files(fastq_dirs):
            usable_count = 0
            i = 0
            print('Correcting barcodes in', transcript_fastq, 'and', barcode_umi_fastq)
            transcript_reader = fastq_reader(transcript_fastq)
            barcode_umi_reader = fastq_reader(barcode_umi_fastq)
            for i, (tr, br) in enumerate(zip(transcript_reader, barcode_umi_reader), 1):
                barcode_pieces = [revcomp(br.seq[s]) for s in reversed(BARCODE_SEGMENTS)]
                corrected = [correcter.correct(barcode) for barcode in barcode_pieces]
                if all(corrected):
                    usable_count += 1
                    umi_seq = br.seq[UMI_SEGMENT]
                    umi_qual = br.qual[UMI_SEGMENT]
                    new_seq = ''.join(corrected + [umi_seq])
                    new_qual = BARCODE_QUAL_DUMMY + umi_qual
                    new_br = Read(
                        read_id=br.read_id,
                        seq=new_seq,
                        unused=br.unused,
                        qual=new_qual,
                    )
                    print(tr.serialize(), file=tro)
                    print(new_br.serialize(), file=cbo)

            print('Total count:', i)
            print('Usable count:', usable_count)
            print('Proportion:', usable_count / i)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('fastq_dirs', type=Path, nargs='+')
    p.add_argument('--barcode_list_file', type=Path, nargs='?')
    args = p.parse_args()

    if args.barcode_list_file is None:
        this_script = Path(__file__)
        args.barcode_list_file = this_script.parent / 'data/snareseq-barcodes.txt'

    main(args.fastq_dirs, args.barcode_list_file)
