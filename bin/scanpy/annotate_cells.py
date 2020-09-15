#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from shutil import copy

from common import Assay

H5AD_PATH = Path('out.h5ad')

def dummy_annotate_cells(h5ad_file: Path):
    copy(h5ad_file, H5AD_PATH)

def main(assay: Assay, h5ad_file: Path):
    if assay == Assay.SCISEQ:
        pass
    elif assay == Assay.SNARESEQ:
        pass
    else:
        print('No annotation to perform for assay', assay)
        dummy_annotate_cells(h5ad_file)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('assay', choices=list(Assay), type=Assay)
    p.add_argument('h5ad_file', type=Path)
    args = p.parse_args()

    main(args.assay, args.h5ad_file)
