#!/usr/bin/env pypy3
from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable

import correct_snareseq_barcodes
import expand_sciseq_barcodes
import extract_slideseq_barcodes

from common import ADJ_OUTPUT_DIR, Assay


def main(assay: Assay, input_dirs: Iterable[Path]):
    ADJ_OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

    if assay == Assay.SCISEQ:
        expand_sciseq_barcodes.main(input_dirs, output_dir=ADJ_OUTPUT_DIR)
    elif assay == Assay.SNARESEQ:
        correct_snareseq_barcodes.main(input_dirs, output_dir=ADJ_OUTPUT_DIR)
    elif assay == Assay.SLIDESEQ:
        extract_slideseq_barcodes.main(input_dirs, output_dir=ADJ_OUTPUT_DIR)
    else:
        print("No barcode adjustment to perform for assay", assay)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("directory", type=Path, nargs="+")
    args = p.parse_args()

    main(args.assay, args.directory)
