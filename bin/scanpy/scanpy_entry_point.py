#!/usr/bin/env python3
from argparse import ArgumentParser
from contextlib import contextmanager
from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

@contextmanager
def new_plot():
    """
    When used in a `with` block, clears matplotlib internal state
    after plotting and saving things. Probably not necessary to be this
    thorough in clearing everything, but extra calls to `plt.clf()` and
    `plf.close()` don't *hurt*

    Intended usage:
        ```
        with new_plot():
            do_matplotlib_things()

            plt.savefig(path)
            # or
            fig.savefig(path)
        ```
    """
    plt.clf()
    try:
        yield
    finally:
        plt.clf()
        plt.close()

def qc_checks(adata: anndata.AnnData):
    qc_by_cell, qc_by_gene = sc.pp.calculate_qc_metrics(adata)

    # current directory is set up by the CWL runner
    qc_path = Path('qc_results.hdf5').absolute()
    print('Saving QC results to', qc_path)
    with pd.HDFStore(qc_path) as store:
        store['qc_by_cell'] = qc_by_cell
        store['qc_by_gene'] = qc_by_gene

def main(h5ad_file: Path):
    adata = anndata.read_h5ad(h5ad_file)
    adata.var_names_make_unique()

    qc_checks(adata)

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis=1)

    adata.raw = sc.pp.log1p(adata, copy=True)

    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

    filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=0.0125, max_mean=5, min_disp=0)
    sc.pl.filter_genes_dispersion(filter_result, show=False, save='dispersion.pdf')

    # %%
    adata = adata[:, filter_result.gene_subset]
    sc.pp.scale(adata, max_value=10)

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
    sc.tl.umap(adata)

    with new_plot():
        sc.pl.umap(adata)
        plt.savefig('umap.pdf', bbox_inches='tight')

    sc.tl.leiden(adata)

    with new_plot():
        sc.pl.umap(adata, color='leiden')
        plt.savefig('umap_by_leiden_cluster.pdf', bbox_inches='tight')

    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

    with new_plot():
        sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
        plt.savefig('marker_genes_by_cluster_t_test.pdf', bbox_inches='tight')

    sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')

    with new_plot():
        sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
        plt.savefig('marker_genes_by_cluster_logreg.pdf', bbox_inches='tight')

    output_file = Path('cluster_marker_genes.h5ad')
    print('Saving output to', output_file.absolute())
    # Save normalized/etc. data
    adata.write_h5ad(output_file)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('alevin_h5ad_file', type=Path)
    args = p.parse_args()

    main(args.alevin_h5ad_file)
