#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import anndata
import matplotlib.pyplot as plt
import scanpy as sc

from common import Assay
from plot_utils import new_plot


def main(assay: Assay, h5ad_file: Path):
    adata = anndata.read_h5ad(h5ad_file)
    if assay.secondary_analysis_layer in adata.layers:
        adata.X = adata.layers[assay.secondary_analysis_layer]
    adata.var_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # add the total counts per cell as observations-annotation to adata
    adata.obs["n_counts"] = adata.X.sum(axis=1)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=5, min_disp=0)

    adata.layers["unscaled"] = adata.X.copy()

    with new_plot():
        sc.pl.highly_variable_genes(adata, show=False)
        plt.savefig("dispersion_plot.pdf", bbox_inches="tight")

    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50)
    sc.tl.umap(adata, min_dist=0.05)
    sc.tl.leiden(adata)

    with new_plot():
        sc.pl.umap(adata, color="leiden", show=False)
        plt.savefig("umap_by_leiden_cluster.pdf", bbox_inches="tight")

    sc.tl.embedding_density(adata, basis="umap")
    with new_plot():
        sc.pl.embedding_density(adata, color_map="viridis_r", show=False)
        plt.savefig("umap_embedding_density.pdf", bbox_inches="tight")

    if "X_spatial" in adata.obsm:
        with new_plot():
            sc.pl.scatter(adata, color="leiden", basis="spatial", show=False)
            plt.savefig("spatial_pos_by_leiden_cluster.pdf", bbox_inches="tight")

    sc.tl.rank_genes_groups(adata, "leiden", method="t-test")

    with new_plot():
        sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
        plt.savefig("marker_genes_by_cluster_t_test.pdf", bbox_inches="tight")

    sc.tl.rank_genes_groups(adata, "leiden", method="logreg")

    with new_plot():
        sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
        plt.savefig("marker_genes_by_cluster_logreg.pdf", bbox_inches="tight")

    output_file = Path("secondary_analysis.h5ad")
    print("Saving output to", output_file.absolute())
    # Save normalized/etc. data
    adata.write_h5ad(output_file)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("alevin_h5ad_file", type=Path)
    args = p.parse_args()

    main(args.assay, args.alevin_h5ad_file)
