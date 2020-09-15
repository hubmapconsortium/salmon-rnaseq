#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from shutil import copy

from common import Assay

H5AD_PATH = Path('out.h5ad')

def dummy_annotate_cells(directory: Path):
    input_h5ad_file = directory / 'out.h5ad'
    copy(input_h5ad_file, H5AD_PATH)

def main(assay: Assay, directory: Path):
    output_path = Path('adj_h5ad')
    output_path.mkdir(exist_ok=True, parents=True)

    if assay == Assay.SCISEQ:
        pass
    elif assay == Assay.SNARESEQ:
        pass
    else:
        dummy_annotate_cells(directory)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('assay', choices=list(Assay), type=Assay)
    p.add_argument('directory', type=Path)
    args = p.parse_args()

    main(args.assay, args.directory)
