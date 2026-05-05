#!/usr/bin/env python3
import re
from argparse import ArgumentParser
from os import fspath, walk
from pathlib import Path
from typing import Iterable, Tuple

import anndata
import manhole
import matplotlib.pyplot as plt
import squidpy as sq
import tifffile as tf

from common import Assay
from plot_utils import new_plot
import spatialdata
import spatialdata_plot

from spatialdata.models import Image2DModel, Image3DModel, Labels2DModel, Labels3DModel, PointsModel, ShapesModel, TableModel
from math import ceil, log2
from aicsimageio import AICSImage
import geopandas
import shapely
import pandas as pd

desired_pixel_size_for_pyramid = 250

ome_tiff_pattern = re.compile(r"(?P<basename>.*)\.ome\.tiff(f?)$")

def find_ome_tiff(input_dir: Path) -> Path:
    for dirpath_str, _, filenames in walk(input_dir):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            if ome_tiff_pattern.match(filename):
                src_filepath = dirpath / filename
                return src_filepath

def get_img_spatialdata(img_dir: Path):
    print("Loading image data")
    image_file = find_ome_tiff(img_dir)
    image = AICSImage(image_file)
    image_data_squeezed = image.data.squeeze()
    print("... done. Original shape:", image.data.shape)
    print(f"New shape: {image_data_squeezed.shape}")

    image_scale_factors = (2,) * ceil(
        log2(max(image_data_squeezed.shape[1:]) / desired_pixel_size_for_pyramid)
    )

    img_for_sdata = Image2DModel.parse(
        data=image_data_squeezed,
        dims=['c', 'y', 'x'],
        scale_factors=image_scale_factors,
    )

    return img_for_sdata


def get_shapes_spatialdata(adata:anndata.AnnData):
    geo_df = geopandas.GeoDataFrame(index=adata.obs.index)
    if 'spatial' in adata.uns:#visium
        radius = adata.uns['spatial']['visium']['scalefactors']['spot_diameter_fullres'] / 2
    else:#slideseq
        radius = 5
        adata.obsm['spatial'] = adata.obsm['X_spatial']
    radius_series = pd.Series(radius, index=geo_df.index)
    coords = adata.obsm['spatial']
    points_list = [shapely.Point(coords[i]) for i in range(len(adata.obs.index))]
    points_series = geopandas.GeoSeries(points_list, index=adata.obs.index)
    geo_df['geometry'] = points_series
    geo_df['radius'] = radius_series
    return ShapesModel.parse(geo_df)

def main(assay: Assay, h5ad_file: Path, img_dir: Path = None):
    if assay in {Assay.VISIUM_FF, Assay.SLIDESEQ}:
        adata = anndata.read_h5ad(h5ad_file)
        # Modify Tissue Coverage Fraction name for squidpy
        adata.obs = adata.obs.rename(columns={'Tissue Coverage Fraction': 'tissue_coverage_fraction'})
        # Instantiate spatialdata_attrs to be able to plot later
        adata.obs['cell_id'] = adata.obs.index
        adata.obs['region'] = 'leiden'
        adata.uns['spatialdata_attrs'] = {'instance_key': 'cell_id', 'region': 'leiden', 'region_key': 'region'}
        # Parse and sanitize anndata object for spatialdata
        adata = spatialdata.sanitize_table(adata, inplace=False)
        table_for_sdata = TableModel.parse(adata)
        # Get shapes
        shapes_for_sdata = get_shapes_spatialdata(adata)
        shapes_for_sdata['leiden'] = table_for_sdata.obs['leiden']
        adata.obsm["spatial"] = adata.obsm["X_spatial"]

        if img_dir:
            tiff_file = find_ome_tiff(input_dir=img_dir)
            img = tf.imread(fspath(tiff_file))
            library_id = list(adata.uns["spatial"].keys())[0]
            adata.uns["spatial"][library_id]["images"] = {"hires": img}
            adata.uns["spatial"][library_id]["scalefactors"] = {
                "tissue_hires_scalef": 1.0,
                "spot_diameter_fullres": 89,
            }

            img_for_sdata = get_img_spatialdata(img_dir)
            sdata = spatialdata.SpatialData(images={'visium_fullres_img':img_for_sdata}, shapes={'visium':shapes_for_sdata}, tables={'table':table_for_sdata})

            sdata.pl.render_images('visium_fullres_img').pl.render_shapes(element='visium', color='leiden').pl.show()
            plt.savefig('spatial_scatter.pdf', bbox_inches='tight')

        else:
            sdata = spatialdata.SpatialData(shapes={'slideseq':shapes_for_sdata}, tables={'table':table_for_sdata})

        output_file_stem_dict = {Assay.VISIUM_FF:"Visium", Assay.SLIDESEQ:"Slideseq"}
        output_file_stem = output_file_stem_dict[assay]
        sdata.write(f'{output_file_stem}.zarr')

        sq.gr.spatial_neighbors(adata)
        sq.gr.nhood_enrichment(adata, cluster_key="leiden")

#        with new_plot():
#            sq.pl.spatial_scatter(adata, color="leiden")
#            plt.savefig("spatial_scatter.pdf", bbox_inches="tight")

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

        #        sq.gr.ripley(adata, cluster_key="leiden")

        #        with new_plot():
        #            sq.pl.ripley(adata, cluster_key="leiden")
        #            plt.savefig("ripley.pdf", bbox_inches="tight")

        output_file = Path("squidpy_annotated.h5ad")
        print("Saving output to", output_file.absolute())
        # Save normalized/etc. data
        # Set column back to expected name
        adata.obs = adata.obs.rename(columns={'tissue_coverage_fraction': 'Tissue Coverage Fraction'})
        adata.write_h5ad(output_file)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("alevin_h5ad_file", type=Path)
    p.add_argument("img_dir", type=Path, nargs="?")

    args = p.parse_args()

    main(args.assay, args.alevin_h5ad_file, args.img_dir)
