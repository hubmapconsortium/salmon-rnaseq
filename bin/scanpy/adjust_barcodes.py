#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import expand_sciseq_barcodes

from common import ADJ_OUTPUT_DIR, Assay

def main(assay: str, input_dir: Path):
    ADJ_OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

    if assay == Assay.SCISEQ:
        expand_sciseq_barcodes.main(input_dir, ADJ_OUTPUT_DIR)
    elif assay == Assay.SNARESEQ:
        pass
    else:
        pass

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('assay', choices=list(Assay), type=Assay)
    p.add_argument('directory', type=Path)
    args = p.parse_args()

    main(args.assay, args.directory)
