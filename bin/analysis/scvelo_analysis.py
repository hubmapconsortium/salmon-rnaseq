#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import anndata
import manhole
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scipy.sparse
import scvelo as scv

from common import Assay
from plot_utils import new_plot

component_neighbor_count = 50
output_file = Path("scvelo_annotated.h5ad")

def main(spliced_h5ad_file: Path, assay: Assay):
    if assay == Assay.VISIUM_FF:
        print("Skipping scVelo analysis for", assay)
        return

    adata = anndata.read_h5ad(spliced_h5ad_file)
    adata.var_names_make_unique()
    print("Before filtering:", adata, sep="\n")

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print("After initial filtering:", adata, sep="\n")

    scv.pp.filter_genes(adata, min_shared_counts=30)
    print("After spliced/unspliced filtering:", adata, sep="\n")
    scv.pp.normalize_per_cell(adata, enforce=True)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
    scv.pp.log1p(adata)

    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50)
    try:
        sc.pp.pca(adata, n_comps=component_neighbor_count)
        sc.pp.neighbors(
            adata,
            n_neighbors=component_neighbor_count,
            n_pcs=component_neighbor_count,
        )
        sc.tl.umap(adata)
        sc.tl.leiden(adata)

        scv.settings.set_figure_params("scvelo")
        scv.utils.show_proportions(adata)

        scv.pp.moments(adata, n_pcs=50, n_neighbors=50)
        scv.pp.moments(
            adata,
            n_pcs=component_neighbor_count,
            n_neighbors=component_neighbor_count,
        )
        scv.tl.recover_dynamics(adata)

        scv.tl.velocity(adata, mode="dynamical")
        scv.tl.velocity_graph(adata)
        scv.tl.velocity_embedding(adata, basis="umap")

        print("Saving output to", output_file)
        adata.write_h5ad(output_file)

        with new_plot():
            scv.pl.velocity_embedding_grid(adata, basis="umap", color="leiden", show=False)
            plt.savefig("scvelo_embedding_grid.pdf", bbox_inches="tight")

    except ValueError as e:
        print("Caught:", e)
        print("Writing empty AnnData file to", output_file)
        empty = anndata.AnnData(X=scipy.sparse.csr_matrix((0, 0), dtype=np.float32))
        empty.write_h5ad(output_file)

        with new_plot():
            message_pieces = [
                "scVelo analysis failed 😵",
                "",
                f"Spliced/unspliced filtering probably kept",
                f"too few cells/nuclei ({adata.shape[0]}) or genes ({adata.shape[1]})",
            ]
            plt.text(0.5, 0.5, "\n".join(message_pieces), horizontalalignment="center")
            plt.savefig("scvelo_embedding_grid.pdf", bbox_inches="tight")


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("alevin_h5ad_file", type=Path)
    p.add_argument("assay", choices=list(Assay), type=Assay)
    args = p.parse_args()

    main(args.alevin_h5ad_file, args.assay)
