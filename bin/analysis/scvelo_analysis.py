#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import anndata
import manhole
import matplotlib.pyplot as plt
import scanpy as sc
import scvelo as scv

from common import Assay
from plot_utils import new_plot


def main(spliced_h5ad_file: Path, assay: Assay):
    if assay not in {Assay.VISIUM_FFPE}:
        adata = anndata.read_h5ad(spliced_h5ad_file)
        adata.var_names_make_unique()

        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)

        scv.pp.filter_genes(adata, min_shared_counts=30)
        scv.pp.normalize_per_cell(adata, enforce=True)

        scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
        scv.pp.log1p(adata)

        sc.pp.pca(adata, n_comps=50)
        sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50)
        sc.tl.umap(adata)
        sc.tl.leiden(adata)

        scv.settings.set_figure_params("scvelo")
        scv.utils.show_proportions(adata)

        scv.pp.moments(adata, n_pcs=50, n_neighbors=50)
        scv.tl.recover_dynamics(adata)

        scv.tl.velocity(adata, mode="dynamical")
        scv.tl.velocity_graph(adata)
        scv.tl.velocity_embedding(adata, basis="umap")

        output_file = Path("scvelo_annotated.h5ad")
        print("Saving output to", output_file)
        adata.write_h5ad(output_file)

        with new_plot():
            scv.pl.velocity_embedding_grid(adata, basis="umap", color="leiden", show=False)
            plt.savefig("scvelo_embedding_grid.pdf", bbox_inches="tight")


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("alevin_h5ad_file", type=Path)
    p.add_argument("assay", choices=list(Assay), type=Assay)
    args = p.parse_args()

    main(args.alevin_h5ad_file, args.assay)
