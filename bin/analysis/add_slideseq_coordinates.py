#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import anndata
import pandas as pd

barcode_matching_dir = "barcode_matching"


def read_bead_pos(dataset_dir: Path) -> pd.DataFrame:
    barcode_matching_dirs = list(dataset_dir.glob(f"**/{barcode_matching_dir}"))
    if l := len(barcode_matching_dirs) != 1:
        message_pieces = [f"Need exactly 1 {barcode_matching_dir} directory, found {l}"]
        message_pieces.extend(f"\t{d}" for d in barcode_matching_dirs)
        raise ValueError("\n".join(message_pieces))

    barcode_dir = barcode_matching_dirs[0]
    all_pos_raw = pd.read_csv(
        barcode_dir / "BeadLocations.txt",
        index_col=None,
        header=None,
    )
    all_pos = all_pos_raw.T

    with open(barcode_dir / "BeadBarcodes.txt") as f:
        all_barcodes = ["".join(line.strip().split(",")) for line in f]
    all_pos.index = all_barcodes
    return all_pos


def annotate(h5ad_path: Path, dataset_dir: Path) -> anndata.AnnData:
    d = anndata.read_h5ad(h5ad_path)
    barcode_pos = read_bead_pos(dataset_dir)

    quant_bc_set = set(d.obs.index)
    pos_bc_set = set(barcode_pos.index)
    overlap = quant_bc_set & pos_bc_set
    positions_overlap = barcode_pos.loc[list(overlap), :]

    quant_minus_pos = quant_bc_set - pos_bc_set
    positions_missing = pd.DataFrame(
        index=list(quant_minus_pos),
        columns=barcode_pos.columns,
        dtype=float,
    )

    quant_pos = pd.concat([positions_overlap, positions_missing])
    quant_pos_ordered = quant_pos.loc[d.obs.index, :]
    d.obsm["X_spatial"] = quant_pos_ordered.to_numpy()

    return d


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("h5ad_path", type=Path)
    p.add_argument("base_dir", type=Path)
    args = p.parse_args()

    d = annotate(args.h5ad_path, args.base_dir)
    d.write_h5ad("slideseq_pos_annotated.h5ad")
