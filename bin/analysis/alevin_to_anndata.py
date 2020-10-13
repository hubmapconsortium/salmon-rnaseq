#!/usr/bin/env python3
from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from anndata import AnnData
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse

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

def sparsify_if_appropriate(mat: scipy.sparse.spmatrix) -> Union[np.ndarray, scipy.sparse.spmatrix]:
    density = mat.nnz / np.prod(mat.shape)
    if density > DENSITY_THRESHOLD:
        return mat.todense()
    else:
        return mat

# noinspection PyDeprecation
def force_dense_matrix(mat: Union[np.ndarray, scipy.sparse.spmatrix]) -> np.matrix:
    # Only used for convenience in doctests
    if isinstance(mat, np.ndarray):
        return np.matrix(mat)
    elif isinstance(mat, scipy.sparse.spmatrix):
        return mat.todense()

@dataclass
class LabeledMatrix:
    matrix: scipy.sparse.spmatrix
    row_labels: List[str]
    col_labels: List[str]

    def to_anndata(self):
        return AnnData(
            X=self.matrix,
            obs=pd.DataFrame(index=self.row_labels),
            var=pd.DataFrame(index=self.col_labels),
        )

def get_col_sum_matrix(
        orig_labels: List[str],
        label_mapping: Dict[str, str]
) -> Tuple[scipy.sparse.spmatrix, List[str]]:
    """
    :param orig_labels:
    :param label_mapping:
    :return: 2-tuple:
      [0] A summation matrix suitable for right-multiplying a data matrix,
          summing columns of that data matrix. Transpose this to sum across rows.
      [1] Labels for the new axis introduced by this summation matrix

    >>> col_labels = list('2143')
    >>> col_mapping = {'2': '3'}
    >>> csm, new_col_labels = get_col_sum_matrix(col_labels, col_mapping)
    >>> csm.todense()
    matrix([[0, 1, 0],
            [1, 0, 0],
            [0, 0, 1],
            [0, 1, 0]])
    >>> new_col_labels
    ['1', '3', '4']
    """
    new_labels = sorted({label_mapping.get(l, l) for l in orig_labels})
    new_label_indices = {label: i for i, label in enumerate(new_labels)}

    data_vec = np.ones(len(orig_labels), dtype=int)
    row_vec = np.arange(len(orig_labels))
    col_vec = np.array([new_label_indices[label_mapping.get(l, l)] for l in orig_labels])

    m = scipy.sparse.coo_matrix((data_vec, (row_vec, col_vec))).tocsr()
    return m, new_labels

def collapse_matrix_rows_cols(
        matrix: LabeledMatrix,
        row_mapping: Optional[Dict[str, str]] = None,
        col_mapping: Optional[Dict[str, str]] = None,
) -> LabeledMatrix:
    """
    Not operating directly on a `AnnData` object, to make it clear that
    this functionality only preserves row and column labels, and does not even
    *attempt* to deal with any supplementary data that could be stored in other
    columns of `AnnData.obs` or `AnnData.var`

    :param matrix: LabeledMatrix instance
    :param row_mapping: Mapping from row labels to new row labels. Any label not
      present in this mapping is returned unchanged.
    :param col_mapping: Mapping from column labels to new column labels. Any label
      not present in this mapping is returned unchanged.
    :return: New LabeledMatrix with entries as sums of appropriate rows and columns

    >>> m = np.arange(1, 10).reshape((3, 3))
    >>> m
    array([[1, 2, 3],
           [4, 5, 6],
           [7, 8, 9]])
    >>> row_labels = list('abc')
    >>> col_labels = list('123')
    >>> lm = LabeledMatrix(matrix=m, row_labels=row_labels, col_labels=col_labels)
    >>> row_mapping = {'a': 'b'}
    >>> col_mapping = {'2': '3'}
    >>> sm = collapse_matrix_rows_cols(lm, row_mapping, col_mapping)
    >>> sm.matrix
    array([[ 5, 16],
           [ 7, 17]])
    >>> sm.row_labels
    ['b', 'c']
    >>> sm.col_labels
    ['1', '3']
    """
    if not (row_mapping or col_mapping):
        return matrix

    if row_mapping is None:
        row_mapping = {}
    if col_mapping is None:
        col_mapping = {}

    col_mat, new_col_labels = get_col_sum_matrix(matrix.col_labels, col_mapping)
    row_mat_t, new_row_labels = get_col_sum_matrix(matrix.row_labels, row_mapping)
    row_mat = row_mat_t.T

    new_data = row_mat @ matrix.matrix @ col_mat
    return LabeledMatrix(
        matrix=new_data,
        row_labels=new_row_labels,
        col_labels=new_col_labels,
    )

def collapse_intron_cols(lm: LabeledMatrix) -> LabeledMatrix:
    intron_mapping = {i: i.split('-')[0] for i in lm.col_labels}
    new_matrix = collapse_matrix_rows_cols(lm, col_mapping=intron_mapping)
    return new_matrix

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

def add_split_spliced_unspliced(lm: LabeledMatrix) -> AnnData:
    """
    :param lm: cell ⨉ gene LabeledMatrix, with columns for both spliced and
      unspliced sequences
    :return: an `AnnData` object (let's say `adata`) constructed from `lm`:
      spliced and unspliced counts are added to compute `adata.X`, spliced
      counts are stored in `adata.layers['spliced']` and unspliced counts
      are likewise stored in `adata.layers['unspliced']`.

    >>> row_labels = ['c1', 'c2']
    >>> col_labels = ['g1', 'g2', 'g2-I', 'g3-I']
    >>> mat = scipy.sparse.csr_matrix(np.arange(1, 9, dtype=np.float32).reshape((2, 4)))
    >>> lm = LabeledMatrix(matrix=mat, row_labels=row_labels, col_labels=col_labels)
    >>> adata: AnnData = add_split_spliced_unspliced(lm)
    >>> list(adata.obs.index)
    ['c1', 'c2']
    >>> list(adata.var.index)
    ['g1', 'g2', 'g3']
    >>> force_dense_matrix(adata.X)
    matrix([[ 1.,  5.,  4.],
            [ 5., 13.,  8.]], dtype=float32)
    >>> force_dense_matrix(adata.layers['spliced'])
    matrix([[1., 2., 0.],
            [5., 6., 0.]], dtype=float32)
    >>> force_dense_matrix(adata.layers['unspliced'])
    matrix([[0., 3., 4.],
            [0., 7., 8.]], dtype=float32)
    """
    # TODO: rethink data types and control flow between helper functions.
    #   These have gone through a few iterations and might need a few cleanups.
    genes_split = [(i, i.split('-')) for i in lm.col_labels]

    all_exons = {g[1][0] for g in genes_split}
    orig_exons = {g[1][0] for g in genes_split if len(g[1]) == 1}
    exons_to_add = list(all_exons - orig_exons)

    d = lm.to_anndata()
    spliced_expanded = expand_anndata(d, list(orig_exons), exons_to_add)

    intron_mapping = {
        g[0]: g[1][0] for g in genes_split
        if len(g[1]) == 2 and g[1][1] == 'I'
    }
    orig_intron_list, mapped_intron_list = (list(z) for z in zip(*intron_mapping.items()))
    introns_to_add = list(all_exons - set(intron_mapping.values()))
    unspliced_expanded = expand_anndata(d, orig_intron_list, introns_to_add, mapped_intron_list)

    collapsed = collapse_intron_cols(lm)

    adata = AnnData(
        X=sparsify_if_appropriate(collapsed.matrix),
        obs=spliced_expanded.obs,
        var=spliced_expanded.var,
        layers={
            'spliced': sparsify_if_appropriate(spliced_expanded.X),
            'unspliced': sparsify_if_appropriate(unspliced_expanded.X),
        },
    )

    return adata

def convert(input_dir: Path) -> Tuple[AnnData, AnnData]:
    """
    :return: 2-tuple:
     [0] full count matrix, with columns for spliced and unspliced regions
     [1] count matrix with `AnnData.X` containing sums of spliced and unspliced,
         `AnnData.layers['spliced']` and `AnnData.layers['unspliced']`
         containing only those counts, respectively
    """
    alevin_dir = input_dir / 'alevin'

    with open(alevin_dir / 'quants_mat_rows.txt') as f:
        cb_names = [line.strip() for line in f]

    with open(alevin_dir / 'quants_mat_cols.txt') as f:
        gene_names = [line.strip() for line in f]

    print('Reading sparse count matrix')
    raw_matrix = scipy.io.mmread(alevin_dir / 'quants_mat.mtx.gz').tocsr()
    raw_labeled = LabeledMatrix(
        matrix=raw_matrix,
        row_labels=cb_names,
        col_labels=gene_names,
    )
    spliced_anndata = add_split_spliced_unspliced(raw_labeled)

    return raw_labeled.to_anndata(), spliced_anndata

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('alevin_output_dir', type=Path)
    args = p.parse_args()

    raw, spliced = convert(args.alevin_output_dir)
    raw.write_h5ad('full.h5ad')
    spliced.write_h5ad('out.h5ad')
