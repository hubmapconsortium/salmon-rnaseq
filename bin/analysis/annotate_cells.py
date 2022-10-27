#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Optional, Sequence

import anndata
import manhole

import add_spatial_coordinates
import annotate_sciseq_barcodes
from common import Assay

H5AD_PATH = Path("expr.h5ad")


def dummy_annotate_cells(h5ad_file: Path) -> anndata.AnnData:
    return anndata.read_h5ad(h5ad_file)


def main(
    assay: Assay,
    h5ad_file: Path,
    raw_fastq_dirs: Sequence[Path],
    metadata_json: Optional[Path],
):
    if assay == Assay.SCISEQ:
        expr_data = annotate_sciseq_barcodes.main(h5ad_file, metadata_json)
    elif assay in {Assay.SLIDESEQ, Assay.VISIUM_FFPE}:
        if len(raw_fastq_dirs) != 1:
            raise ValueError("Need exactly 1 input directory for Slide-seq")
        expr_data = add_spatial_coordinates.annotate(h5ad_file, raw_fastq_dirs[0], assay)
    else:
        print("No annotation to perform for assay", assay)
        expr_data = dummy_annotate_cells(h5ad_file)

    expr_data.write_h5ad(H5AD_PATH)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("h5ad_file", type=Path)
    p.add_argument("raw_fastq_dir", type=Path, nargs="+")
    p.add_argument("--metadata_json", type=Path)
    args = p.parse_args()

    main(args.assay, args.h5ad_file, args.raw_fastq_dir, args.metadata_json)
