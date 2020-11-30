#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Optional

import anndata
import annotate_sciseq_barcodes

from common import Assay

H5AD_PATH = Path("expr.h5ad")
ZARR_PATH = H5AD_PATH.with_suffix(".zarr")


def dummy_annotate_cells(h5ad_file: Path) -> anndata.AnnData:
    return anndata.read_h5ad(h5ad_file)


def main(assay: Assay, h5ad_file: Path, metadata_json: Optional[Path]):
    if assay == Assay.SCISEQ:
        expr_data = annotate_sciseq_barcodes.main(h5ad_file, metadata_json)
    else:
        print("No annotation to perform for assay", assay)
        expr_data = dummy_annotate_cells(h5ad_file)

    expr_data.write_h5ad(H5AD_PATH)
    expr_data.write_zarr(ZARR_PATH)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("h5ad_file", type=Path)
    p.add_argument("--metadata_json", type=Path)
    args = p.parse_args()

    main(args.assay, args.h5ad_file, args.metadata_json)
