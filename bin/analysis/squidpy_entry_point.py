#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from os import walk, fspath
import anndata
import manhole
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

from common import Assay
from plot_utils import new_plot
from typing import Iterable, Tuple
import re

import cv2
import numpy as np

ome_tiff_pattern = re.compile(r"(?P<basename>.*)\.ome\.tiff(f?)$")

def find_ome_tiffs(input_dir: Path) -> Iterable[Path]:
    """
    Yields 2-tuples:
     [0] full Path to source file
     [1] output file Path (source file relative to input_dir)
    """
    for dirpath_str, _, filenames in walk(input_dir):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            if ome_tiff_pattern.match(filename):
                src_filepath = dirpath / filename
                yield src_filepath

def main(assay: Assay, h5ad_file: Path, img_dir: Path = None):
    if assay in {Assay.VISIUM_FF, Assay.SLIDESEQ}:
        adata = anndata.read(h5ad_file)
        adata.obsm["spatial"] = adata.obsm["X_spatial"]
        if img_dir:
            tiff_file = list(find_ome_tiffs(input_dir=img_dir))[0]
#            img = cv2.imread(fspath(tiff_file[0]))

        sq.gr.spatial_neighbors(adata)
        sq.gr.nhood_enrichment(adata, cluster_key="leiden")

        with new_plot():
            sq.pl.nhood_enrichment(adata, cluster_key="leiden")
            plt.savefig("neighborhood_enrichment.pdf", bbox_inches="tight")

        sq.gr.co_occurrence(adata, cluster_key="leiden")

        with new_plot():
            sq.pl.co_occurrence(adata, cluster_key="leiden")
            plt.savefig("co_occurrence.pdf", bbox_inches="tight")

        sq.gr.centrality_scores(adata, cluster_key="leiden")

        with new_plot():
            sq.pl.centrality_scores(adata, cluster_key="leiden")
            plt.savefig("centrality_scores.pdf", bbox_inches="tight")

        sq.gr.interaction_matrix(adata, cluster_key="leiden")

        with new_plot():
            sq.pl.interaction_matrix(adata, cluster_key="leiden")
            plt.savefig("interaction_matrix.pdf", bbox_inches="tight")

        sq.gr.ripley(adata, cluster_key="leiden")

        with new_plot():
            sq.pl.ripley(adata, cluster_key="leiden")
            plt.savefig("ripley.pdf", bbox_inches="tight")

        output_file = Path("squidpy_annotated.h5ad")
        print("Saving output to", output_file.absolute())
        # Save normalized/etc. data
        adata.write_h5ad(output_file)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("alevin_h5ad_file", type=Path)
    p.add_argument("img_dir", type=Path, nargs='?')

    args = p.parse_args()

    main(args.assay, args.alevin_h5ad_file, args.img_dir)
