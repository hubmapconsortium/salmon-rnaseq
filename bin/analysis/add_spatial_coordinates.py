#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import anndata
import manhole
import pandas as pd
from common import Assay
from os import walk
from typing import Iterable

barcode_matching_dir = "barcode_matching"

def read_slideseq_pos(dataset_dir: Path) -> pd.DataFrame:
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


def find_files(directory: Path, pattern: str) -> Iterable[Path]:
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            if filepath.match(pattern):
                yield filepath

def read_visium_pos(dataset_dir: Path) -> pd.DataFrame:
    gpr_file = list(find_files(dataset_dir, "*.gpr"))[0]
    gpr_df = pd.read_csv(gpr_file, sep='\t', skiprows=9)
    print(f"len gpr_df.index: {len(gpr_df.index)}")
    gpr_df = gpr_df[gpr_df['Block'] == 1]
    gpr_df = gpr_df.set_index(['Column', 'Row'], inplace=False, drop=True)
    plate_version_number = gpr_file.stem[1]
    barcode_coords_file = Path(f"/opt/data/visium-v{plate_version_number}_coordinates.txt")
    coords_df = pd.read_csv(barcode_coords_file, sep='\t', names=['barcode', 'Column', 'Row'])
    coords_df = coords_df.set_index(['Column', 'Row'])
    print(f"len coords_df.index: {len(coords_df.index)}")
    gpr_df['barcode'] = coords_df['barcode']
    gpr_df = gpr_df[['barcode', 'X', 'Y']]
    print(f"len gpr_df with barcodes: {len(gpr_df.index)}")
    gpr_df = gpr_df[~gpr_df.barcode.isna()]
    print(f"len gpr_df with barcodes non_null: {len(gpr_df.index)}")
    gpr_df = gpr_df.reset_index(inplace=False)
    gpr_df = gpr_df.set_index('barcode', inplace=False, drop=True)
    return gpr_df

def annotate(h5ad_path: Path, dataset_dir: Path, assay: Assay) -> anndata.AnnData:
    assert assay in {Assay.SLIDESEQ, Assay.VISIUM_FFPE}
    d = anndata.read_h5ad(h5ad_path)
    print(f"adata.obs.index: {len(d.obs.index)}")
    if assay == Assay.SLIDESEQ:
        barcode_pos = read_slideseq_pos(dataset_dir)
    elif assay in {Assay.VISIUM_FFPE, Assay.VISIUM_FF}:
        barcode_pos = read_visium_pos(dataset_dir)

    quant_bc_set = set(d.obs.index)
    print(f"len quant_bc set: {len(quant_bc_set)}")
    pos_bc_set = set(barcode_pos.index)
    print(f"len pos_bc set:  {len(pos_bc_set)}")
    overlap = quant_bc_set & pos_bc_set
    print(f"len overlap: {len(overlap)}")
    positions_overlap = barcode_pos.loc[list(overlap), :]

    quant_minus_pos = quant_bc_set - pos_bc_set
    print(f"len quant_minus_pos: {len(quant_minus_pos)}")
    print(f"len pos_minus_quant: {len(pos_bc_set - quant_bc_set)}")
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
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("h5ad_path", type=Path)
    p.add_argument("base_dir", type=Path)
    p.add_argument("assay", choices=list(Assay), type=Assay)
    args = p.parse_args()

    d = annotate(args.h5ad_path, args.base_dir, args.assay)
    d.write_h5ad("spatial_pos_annotated.h5ad")
