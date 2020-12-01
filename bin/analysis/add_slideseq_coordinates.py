#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import anndata
import pandas as pd


def read_bead_pos(fastq_dir: Path) -> pd.DataFrame:
    all_pos_raw = pd.read_csv(
        fastq_dir / "BeadLocations.txt",
        index_col=None,
        header=None,
    )
    all_pos = all_pos_raw.T

    with open(fastq_dir / "BeadBarcodes.txt") as f:
        all_barcodes = ["".join(line.strip().split(",")) for line in f]
    all_pos.index = all_barcodes
    return all_pos


def annotate(h5ad_path: Path, raw_fastq_dir: Path) -> anndata.AnnData:
    d = anndata.read_h5ad(h5ad_path)
    barcode_pos = read_bead_pos(raw_fastq_dir)

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
    d.obsm["X_slideseq"] = quant_pos_ordered.to_numpy()

    return d


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("h5ad_path", type=Path)
    p.add_argument("base_dir", type=Path)
    args = p.parse_args()

    d = annotate(args.h5ad_path, args.base_dir)
    d.write_h5ad("slideseq_pos_annotated.h5ad")
