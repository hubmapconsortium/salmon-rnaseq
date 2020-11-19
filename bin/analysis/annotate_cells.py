#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from shutil import copy
from typing import Optional

import annotate_sciseq_barcodes

from common import Assay

H5AD_PATH = Path("expr.h5ad")


def dummy_annotate_cells(h5ad_file: Path):
    copy(h5ad_file, H5AD_PATH)


def main(assay: Assay, h5ad_file: Path, metadata_json: Optional[Path]):
    if assay == Assay.SCISEQ:
        annotate_sciseq_barcodes.main(h5ad_file, metadata_json)
    else:
        print("No annotation to perform for assay", assay)
        dummy_annotate_cells(h5ad_file)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("h5ad_file", type=Path)
    p.add_argument("--metadata_json", type=Path)
    args = p.parse_args()

    main(args.assay, args.h5ad_file, args.metadata_json)
