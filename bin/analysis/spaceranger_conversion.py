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
    filtered_file = spaceranger_dir / Path("filtered_feature_bc_matrix.h5")
    raw_file = spaceranger_dir / Path("raw_feature_bc_matrix.h5")

    if assay in {Assay.VISIUM_FF}:
        filtered_adata = sq.read_visium(filtered_file)
        filtered_adata.write("filtered_spaceranger.h5ad")
        raw_adata = sq.read_visium(raw_file)
        raw_adata.write("raw_spaceranger.h5ad")

if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("spaceranger_dir", type=Path, nargs='?')

    args = p.parse_args()

    main(args.assay, args.alevin_h5ad_file, args.spaceranger_dir)
