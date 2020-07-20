#!/usr/bin/env python3
from argparse import ArgumentParser
from functools import lru_cache
from pathlib import Path
import re
from shutil import copy
from typing import Dict

from utils import Read, fastq_reader, smart_open

BARCODE_DATA_DIR = here = Path(__file__).parent / 'data'
FASTQ_INPUT_PATTERN = re.compile(r'(?P<basename>.+)\.(fastq|fq)(.gz)?')

BASE_OUTPUT_DIR = Path('adj_fastq')

@lru_cache(maxsize=None)
def get_base_qual_str(length: int) -> str:
    return 'F' * length

class BarcodeMapper:
    p5_mapping: Dict[str, str]
    p7_mapping: Dict[str, str]
    rt2_mapping: Dict[str, str]
    rt2_mapping: Dict[str, str]

    @classmethod
    def read_barcode_mapping(cls, filename) -> Dict[str, str]:
        mapping = {}
        with open(BARCODE_DATA_DIR / filename) as f:
            for line in f:
                pieces = line.strip().split()
                mapping[pieces[0]] = pieces[1]
        return mapping

    def __init__(self):
        labels = ['p5', 'p7', 'rt', 'rt2']
        for label in labels:
            setattr(self, f'{label}_mapping', self.read_barcode_mapping(f'{label}.txt'))

def convert(mapper: BarcodeMapper, input_fastq: Path, basename: str):
    print('Converting', input_fastq)
    barcode_umi_path = BASE_OUTPUT_DIR / f'{basename}_R1.fastq'
    transcript_path = BASE_OUTPUT_DIR / f'{basename}_R2.fastq.gz'
    copy(input_fastq, transcript_path)

    with smart_open(barcode_umi_path, 'wt') as f:
        for transcript_read in fastq_reader(input_fastq):
            id_pieces = transcript_read.read_id.lstrip('@').split('|')
            p7 = mapper.p7_mapping[id_pieces[2]]
            p5 = mapper.p5_mapping[id_pieces[3]]
            rt2 = mapper.rt2_mapping[id_pieces[4]]
            umi = id_pieces[5]

            barcode_umi = p7 + p5 + rt2 + umi
            barcode_umi_read = Read(
                read_id=transcript_read.read_id,
                seq=barcode_umi,
                unused=transcript_read.unused,
                qual=get_base_qual_str(len(barcode_umi)),
            )
            print(barcode_umi_read.serialize(), file=f)

def main(directory: Path):
    mapper = BarcodeMapper()
    for child in directory.iterdir():
        if m := FASTQ_INPUT_PATTERN.match(child.name):
            convert(mapper, child, m.group('basename'))

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('directory', type=Path)
    args = p.parse_args()

    main(args.directory)
