#!/usr/bin/env python3
import json
from argparse import ArgumentParser
from pathlib import Path
from shutil import copy
from typing import Dict, List, Optional, Sequence, Tuple, Union

import manhole
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
from anndata import AnnData
from fastq_utils import smart_open

from common import AnnDataLayer, Assay

DATA_PATH = Path("/opt/data")
DEFAULT_HUGO_ENSEMBL_MAPPING_PATH = DATA_PATH / "ensembl_hugo_mapping.json.xz"
GENOME_BUILD_PATH = DATA_PATH / "genome_build.json"

# As per https://gist.github.com/flying-sheep/f46e89b388fed736ff0b68fb8fd83af6
# the break-even point for density seems to be around 0.6 to 0.7 for large enough
# data sets (and ours are definitely large enough). Past that point, sparse
# formats actually require *more* storage.
#
# Even though density between 0.5 and ~0.65 still shows space savings with a
# sparse format, there is computational overhead in working with sparse matrices,
# so let's set our cutoff a little below 0.65.
#
# If the density of the data set is <= this, store as a SciPy CSR sparse matrix:
DENSITY_THRESHOLD = 0.5


def sparsify_if_appropriate(
    mat: scipy.sparse.spmatrix,
) -> Union[np.ndarray, scipy.sparse.spmatrix]:
    density = mat.nnz / np.prod(mat.shape)
    if density > DENSITY_THRESHOLD:
        return mat.todense()
    else:
        return mat


class EnsemblHugoMapper:
    mapping: Dict[str, str]

    def __init__(self):
        self.mapping = {}

    @classmethod
    def populate(cls, path: Path, assay: Assay):
        self = cls()
        with smart_open(path) as f:
            self.mapping.update(json.load(f))
        return self

    def annotate(self, data: AnnData, assay: Assay):

        symbols = [self.mapping.get(e) for e in data.var.index]
        data.var.loc[
            :,
            "hugo_symbol",
        ] = symbols


def build_anndata(X, rows: Sequence[str], cols: Sequence[str], **kwargs) -> AnnData:
    """
    Helper to construct an AnnData object from a data matrix and
    specified rows/columns, with calls to pd.DataFrame for obs and
    var data structures
    """
    return AnnData(
        X=X,
        obs=pd.DataFrame(index=rows),
        var=pd.DataFrame(index=cols),
        **kwargs,
    )


# noinspection PyDeprecation
def force_dense_matrix(mat: Union[np.ndarray, scipy.sparse.spmatrix]) -> np.matrix:
    # Only used for convenience in doctests
    if isinstance(mat, np.ndarray):
        return np.matrix(mat)
    elif isinstance(mat, scipy.sparse.spmatrix):
        return mat.todense()


def expand_anndata(
    d: AnnData,
    selected_cols: List[str],
    added_cols: List[str],
    replacement_selected_cols: Optional[List] = None,
) -> AnnData:
    x_real = d[:, selected_cols].X
    x_addition = scipy.sparse.coo_matrix(
        (d.shape[0], len(added_cols)),
        dtype=d.X.dtype,
    )
    x_full = scipy.sparse.hstack([x_real, x_addition])
    orig_cols = replacement_selected_cols or selected_cols
    expanded_cols = orig_cols + added_cols
    d_full = AnnData(
        X=x_full.tocsr(),
        obs=d.obs.copy(),
        var=pd.DataFrame(index=expanded_cols),
    )
    d_full_sorted = d_full[:, sorted(expanded_cols)].copy()
    return d_full_sorted


def add_split_spliced_unspliced(d: AnnData) -> AnnData:
    """
    :param d: AnnData with columns for both spliced and unspliced sequences
    :return: an `AnnData` object (let's say `adata`) constructed from `d`:
      spliced counts are stored in `adata.X` and `adata.layers['spliced']`,
      unspliced counts are stored in `adata.layers['unspliced']`, and the
      sum is stored in `adata.layers['spliced_unspliced_sum']`.

    >>> row_labels = ['c2', 'c1']
    >>> col_labels = ['g2', 'g1', 'g3-I', 'g2-I']
    >>> mat = scipy.sparse.csr_matrix(np.arange(1, 9, dtype=np.float32).reshape((2, 4)))
    >>> d = build_anndata(X=mat, rows=row_labels, cols=col_labels)
    >>> adata: AnnData = add_split_spliced_unspliced(d)
    >>> list(adata.obs.index)
    ['c1', 'c2']
    >>> list(adata.var.index)
    ['g1', 'g2', 'g3']
    >>> force_dense_matrix(adata.X)
    matrix([[6., 5., 0.],
            [2., 1., 0.]], dtype=float32)
    >>> force_dense_matrix(adata.layers['spliced'])
    matrix([[6., 5., 0.],
            [2., 1., 0.]], dtype=float32)
    >>> force_dense_matrix(adata.layers['unspliced'])
    matrix([[0., 8., 7.],
            [0., 4., 3.]], dtype=float32)
    >>> force_dense_matrix(adata.layers['spliced_unspliced_sum'])
    matrix([[ 6., 13.,  7.],
            [ 2.,  5.,  3.]], dtype=float32)
    """
    d_sorted = d[sorted(d.obs.index), sorted(d.var.index)]
    # TODO: rethink data types and control flow between helper functions.
    #   These have gone through a few iterations and might need a few cleanups.
    genes_split = [(i, i.split("-")) for i in d_sorted.var.index]

    all_exons = {g[1][0] for g in genes_split}
    orig_exons = {g[1][0] for g in genes_split if len(g[1]) == 1}
    exons_to_add = list(all_exons - orig_exons)

    spliced_expanded = expand_anndata(d_sorted, list(orig_exons), exons_to_add)

    intron_mapping = {g[0]: g[1][0] for g in genes_split if len(g[1]) == 2 and g[1][1] == "I"}
    orig_intron_list, mapped_intron_list = (list(z) for z in zip(*intron_mapping.items()))
    introns_to_add = list(all_exons - set(intron_mapping.values()))
    unspliced_expanded = expand_anndata(
        d_sorted,
        orig_intron_list,
        introns_to_add,
        mapped_intron_list,
    )

    spliced = sparsify_if_appropriate(spliced_expanded.X)
    adata = AnnData(
        X=spliced,
        obs=spliced_expanded.obs,
        var=spliced_expanded.var,
        layers={
            AnnDataLayer.SPLICED: spliced,
            AnnDataLayer.UNSPLICED: sparsify_if_appropriate(unspliced_expanded.X),
            AnnDataLayer.SPLICED_UNSPLICED_SUM: sparsify_if_appropriate(
                spliced_expanded.X + unspliced_expanded.X
            ),
        },
    )

    return adata


def convert(
    assay: Assay, input_dir: Path, ensembl_hugo_mapping_path: Path, organism: Optional[str]="human",
) -> Tuple[AnnData, AnnData]:
    """
    :return: 2-tuple:
     [0] full count matrix, with columns for spliced and unspliced regions
     [1] count matrix with `AnnData.X` containing sums of spliced and unspliced,
         `AnnData.layers['spliced']` and `AnnData.layers['unspliced']`
         containing only those counts, respectively
    """
    alevin_dir = input_dir / "alevin"

    with open(alevin_dir / "quants_mat_rows.txt") as f:
        cb_names = [line.strip() for line in f]

    with open(alevin_dir / "quants_mat_cols.txt") as f:
        gene_names = [line.strip() for line in f]

    print("Reading sparse count matrix")
    raw_matrix = scipy.io.mmread(alevin_dir / "quants_mat.mtx.gz").tocsr()
    raw_labeled = build_anndata(
        X=raw_matrix,
        rows=cb_names,
        cols=gene_names,
    )

    spliced_anndata = add_split_spliced_unspliced(raw_labeled)

    ensembl_hugo_mapper = EnsemblHugoMapper.populate(ensembl_hugo_mapping_path, assay)
    ensembl_hugo_mapper.annotate(raw_labeled, assay)
    ensembl_hugo_mapper.annotate(spliced_anndata, assay)

    return raw_labeled, spliced_anndata


if __name__ == "__main__":
    manhole.install(activate_on="USR1")
    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("alevin_output_dir", type=Path)
    p.add_argument(
        "--ensembl_hugo_mapping_path",
        type=Path,
        default=DEFAULT_HUGO_ENSEMBL_MAPPING_PATH,
    )
    p.add_argument("--organism", type=str, nargs="?", default="human")

    args = p.parse_args()

    raw, spliced = convert(args.assay, args.alevin_output_dir, args.ensembl_hugo_mapping_path, args.organism)
    if raw:
        raw.write_h5ad("raw_expr.h5ad")
    print(spliced)
    spliced.write_h5ad("expr.h5ad")
    if GENOME_BUILD_PATH.is_file():
        copy(GENOME_BUILD_PATH, Path.cwd() / GENOME_BUILD_PATH.name)
