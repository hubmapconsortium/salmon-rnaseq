#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Dict, Optional

import anndata
import pandas as pd

BARCODE_DATA_DIR = here = Path(__file__).parent / 'data'

BARCODE_LENGTH = 10
BARCODE_STARTS = [0, 10, 20]
BARCODE_SEGMENTS = [slice(start, start + BARCODE_LENGTH) for start in BARCODE_STARTS]
UMI_SEGMENT = slice(0, 8)

class CellIdMapper:
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
                mapping[pieces[1]] = pieces[0]
        return mapping

    def __init__(self):
        labels = ['p5', 'p7', 'rt', 'rt2']
        for label in labels:
            setattr(self, f'{label}_mapping', self.read_barcode_mapping(f'{label}.txt'))

def annotate(mapper: CellIdMapper, h5ad_file: Path, fastq_basename: Optional[str] = None):
    d = anndata.read_h5ad(h5ad_file)
    rt2s = []
    cell_ids = []
    for i in d.obs.index:
        barcode_pieces = [i[seg] for seg in BARCODE_SEGMENTS]
        p7 = mapper.p7_mapping[barcode_pieces[0]]
        p5 = mapper.p5_mapping[barcode_pieces[1]]
        rt2s.append(mapper.rt2_mapping[barcode_pieces[2]])

        sample_basename = f'{fastq_basename}-' if fastq_basename else ''
        cell_id = f'{sample_basename}P7{p7}-P5{p5}'
        cell_ids.append(cell_id)

    d.obs.loc[:, 'cell_id'] = pd.Series(cell_ids, index=d.obs.index).astype('string')
    d.obs.loc[:, 'sample'] = pd.Series(rt2s, index=d.obs.index).astype('category')

    d.write_h5ad('out.h5ad')

def main(h5ad_file: Path, fastq_basename: Optional[str]):
    mapper = CellIdMapper()
    annotate(mapper, h5ad_file, fastq_basename)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('h5ad_file', type=Path)
    p.add_argument('--fastq_basename')
    args = p.parse_args()

    main(args.h5ad_file, args.fastq_basename)
