#!/usr/bin/env python3
from functools import lru_cache
from pathlib import Path
import re
from shutil import copy
from typing import Dict

from fastq_utils import Read, fastq_reader, smart_open

# Relative to path *in container*
BARCODE_DATA_DIR = Path(__file__).parent / 'data/sciseq'
FASTQ_INPUT_PATTERN = re.compile(r'(?P<basename>.+)\.(fastq|fq)(.gz)?')

BASE_OUTPUT_DIR = Path('adj_fastq')

@lru_cache(maxsize=None)
def get_base_qual_str(length: int) -> str:
    return 'F' * length

class BarcodeMapper:
    p5_mapping: Dict[str, str]
    p7_mapping: Dict[str, str]
    rt_mapping: Dict[str, str]
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

def convert(mapper: BarcodeMapper, input_fastq: Path, output_dir: Path, basename: str):
    output_dir.mkdir(exist_ok=True, parents=True)
    print('Converting', input_fastq)
    barcode_umi_path = output_dir / f'{basename}_R1.fastq'
    transcript_path = output_dir / f'{basename}_R2.fastq.gz'
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

def main(directory: Path, output_dir):
    mapper = BarcodeMapper()
    for child in directory.iterdir():
        if m := FASTQ_INPUT_PATTERN.match(child.name):
            convert(mapper, child, output_dir, m.group('basename'))
