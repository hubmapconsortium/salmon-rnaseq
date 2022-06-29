#!/usr/bin/env python3
import json
from argparse import ArgumentParser
from pathlib import Path
from typing import List, Optional, Tuple

import anndata
import manhole
import pandas as pd
import scanpy as sc

from common import Assay

alevin_keys_to_write: List[Tuple[str, Optional[str]]] = [
    ("total_reads", None),
    ("noisy_cb_reads", None),
    ("noisy_umi_reads", None),
    ("used_reads", None),
    ("mapping_rate", "pct_reads_mapped"),
    ("total_cbs", "barcodes_total"),
    ("used_cbs", "barcodes_used"),
    ("final_num_cbs", "barcodes_pre_filtering"),
    ("deduplicated_umis", None),
    ("mean_umis_per_cell", None),
    ("mean_genes_per_cell", None),
]


def read_salmon_output(salmon_dir: Path):
    aux_info_path = salmon_dir / "aux_info/alevin_meta_info.json"
    with open(aux_info_path) as f:
        aux_info = json.load(f)

    return aux_info


def write_scanpy_qc(adata: anndata.AnnData):
    qc_by_cell, qc_by_gene = sc.pp.calculate_qc_metrics(adata)

    qc_path = Path("qc_results.hdf5").absolute()
    print("Saving QC results to", qc_path)
    with pd.HDFStore(qc_path) as store:
        store["qc_by_cell"] = qc_by_cell
        store["qc_by_gene"] = qc_by_gene


def main(assay: Assay, h5ad_primary: Path, h5ad_secondary: Path, salmon_dir: Path):
    expr_primary = anndata.read_h5ad(h5ad_primary)
    if assay.secondary_analysis_layer in expr_primary.layers:
        expr_primary.X = expr_primary.layers[assay.secondary_analysis_layer]
    expr_primary.var_names_make_unique()

    write_scanpy_qc(expr_primary)

    spliced_total = expr_primary.layers["spliced"].sum()
    unspliced_total = expr_primary.layers["unspliced"].sum()
    spliced_unspliced_total = expr_primary.layers["spliced_unspliced_sum"].sum()

    expr_secondary = anndata.read_h5ad(h5ad_secondary)

    salmon_measures = read_salmon_output(salmon_dir)

    data_to_write = {
        (maybe_new_name or key): salmon_measures[key]
        for key, maybe_new_name in alevin_keys_to_write
    }

    data_to_write["barcodes_pre_filtering"] = expr_secondary.shape[0]
    data_to_write["pct_counts_exonic"] = (spliced_total / spliced_unspliced_total) * 100
    data_to_write["pct_counts_intronic"] = (unspliced_total / spliced_unspliced_total) * 100

    with open("qc_results.json", "w") as f:
        json.dump(data_to_write, f, indent=4)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")
    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("h5ad_primary", type=Path)
    p.add_argument("h5ad_secondary", type=Path)
    p.add_argument("salmon_dir", type=Path)
    args = p.parse_args()

    main(args.assay, args.h5ad_primary, args.h5ad_secondary, args.salmon_dir)
