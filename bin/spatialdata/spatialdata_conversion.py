#!/usr/bin/env python3
import re
from argparse import ArgumentParser
from os import fspath, walk
from pathlib import Path
from typing import Iterable, Tuple

import anndata
from geopandas import GeoSeries, GeoDataFrame
from shapely.geometry import Point
import tifffile as tf
import manhole
from spatialdata import SpatialData

from common import Assay

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

def get_cell_circles(adata, diameter):
    a = adata.obsm["X_spatial"]
    points = [Point(a[[i],:]) for i in range(len(adata.obs.index))]
    points_series = GeoSeries(data=points, index=adata.obs.index)
    diameters_series = GeoSeries(data=[diameter for i in adata.obs.index], index=adata.obs.index)
    geo_df = GeoDataFrame(index=adata.obs.index)
    geo_df['geometry'] = points_series
    geo_df['diameter'] = diameters_series
    return geo_df

def main(assay: Assay, h5ad_file: Path, img_dir: Path = None):
    if assay in {Assay.VISIUM_FF}:
        adata = anndata.read(h5ad_file)
        adata.obsm["spatial"] = adata.obsm["X_spatial"]
        if img_dir:
            tiff_file = list(find_ome_tiffs(input_dir=img_dir))[0]
            img = tf.imread(fspath(tiff_file))
            library_id = list(adata.uns["spatial"].keys())[0]
            spot_diameter = adata.uns["spatial"][library_id]["spot_diameter_fullres"]
            cell_circles = get_cell_circles(adata, spot_diameter)
            sdata = SpatialData()
            sdata['fullres'] = img
            sdata['capture_beads'] = cell_circles
            sdata['table'] = adata
            sdata.write('Visium.zarr')


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("alevin_h5ad_file", type=Path)
    p.add_argument("img_dir", type=Path, nargs="?")

    args = p.parse_args()

    main(args.assay, args.alevin_h5ad_file, args.img_dir)
