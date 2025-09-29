#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import anndata as ad
import pandas as pd

from DeepScence.api import DeepScence


def main(h5ad_file):
    adata = ad.read_h5ad(h5ad_file)

    # DeepScence is expecting HUGO symbols in the var index
    adata.var['ensembl_id'] = adata.var.index
    adata.var['new_index'] = adata.var['hugo_symbol'].astype(str)
    adata.var['new_index'] = adata.var['new_index'].fillna(adata.var['ensembl_id'])
    adata.var.index = adata.var['new_index']

    # Run DeepScence
    adata = DeepScence(adata, binarize=True)

    # Rename and reorganize new columns
    df = pd.DataFrame()
    df['ds'] = adata.obs['ds']
    df['binary'] = adata.obs['binary']
    df['SnC'] = adata.obs['SnC']
    adata.obsm['DeepScence'] = df
    adata.obs = adata.obs.drop(['ds', 'binary', 'SnC'], axis=1)
    adata.var.index = adata.var['ensembl_id']
    adata.var = adata.var.drop(['new_index', 'ensembl_id'], axis=1)
    adata.uns['DeepScence'] = "Qu Y, Ji B, Dong R, et al. Single-cell and spatial detection of senescent cells using DeepScence. Preprint. bioRxiv. 2025;2023.11.21.568150. Published 2025 Apr 23. doi:10.1101/2023.11.21.568150"
    adata.uns['DeepScence_log'] = adata.uns['log']
    del adata.uns['log']

    adata.write_h5ad("expr.h5ad")


if __name__ == "__main__":

    p = ArgumentParser()
    p.add_argument("h5ad_file", type=Path)
    args = p.parse_args()

    main(
        args.h5ad_file
    )