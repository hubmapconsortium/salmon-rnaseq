#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

from os import walk, fspath
import anndata
import manhole
from pathlib import Path
import squidpy as sq
from common import Assay


def main(assay: Assay, spaceranger_dir: Path = None):

    if assay in {Assay.VISIUM_FF}:
        filtered_adata = sq.read.visium(spaceranger_dir, counts_file="filtered_feature_bc_matrix.h5")
        filtered_adata.write("filtered_spaceranger.h5ad")
        raw_adata = sq.read.visium(spaceranger_dir, counts_file="raw_feature_bc_matrix.h5")
        raw_adata.write("raw_spaceranger.h5ad")

if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("spaceranger_dir", type=Path, nargs='?')

    args = p.parse_args()

    main(args.assay, args.spaceranger_dir)
