#!/usr/bin/env python3
from functools import lru_cache
import json
from pathlib import Path
import re
from subprocess import run
from typing import Dict

from fastq_utils import Read, fastq_reader, smart_open

# Relative to path *in container*
BARCODE_DATA_DIR = Path(__file__).parent / 'data/sciseq'
FASTQ_INPUT_PATTERN = re.compile(r'(?P<basename>.+)\.(fastq|fq)(.gz)?')

BASE_OUTPUT_DIR = Path('adj_fastq')
METADATA_JSON_PATH = Path('metadata.json')

# Probably faster than piping through the Python interpreter, even
# though we're reading everything anyway, to write the barcode/UMI FASTQ
GUNZIP_COMMAND = ['gunzip', '-c', '{path}']

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

def decompress_fastq(input_path: Path, output_path: Path):
    print('Decompressing', input_path, 'to', output_path)
    with open(output_path, 'wb') as o:
        command = [
            piece.format(path=input_path)
            for piece in GUNZIP_COMMAND
        ]
        run(command, stdout=o, check=True)
    return output_path

def convert(mapper: BarcodeMapper, input_fastq: Path, output_dir: Path, basename: str):
    output_dir.mkdir(exist_ok=True, parents=True)
    print('Converting', input_fastq)
    barcode_umi_path = output_dir / f'{basename}_R1.fastq'
    transcript_path = output_dir / f'{basename}_R2.fastq'
    decompress_fastq(input_fastq, transcript_path)

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
    experiment_ids = set()
    for child in directory.iterdir():
        if m := FASTQ_INPUT_PATTERN.match(child.name):
            basename = m.group('basename')
            convert(mapper, child, output_dir, basename)
            experiment_ids.add(basename.split('-')[0])

    # TODO: relax this
    assert len(experiment_ids) == 1
    experiment_id = next(iter(experiment_ids))
    with open(METADATA_JSON_PATH, 'w') as f:
        json.dump({'experiment_id': experiment_id}, f)
