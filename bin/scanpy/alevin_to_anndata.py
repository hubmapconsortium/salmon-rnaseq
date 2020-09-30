#!/usr/bin/env python3
from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata
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

@dataclass
class LabeledMatrix:
    matrix: scipy.sparse.spmatrix
    row_labels: List[str]
    col_labels: List[str]

    def to_anndata(self):
        return anndata.AnnData(
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
    Not operating directly on a `anndata.AnnData` object, to make it clear that
    this functionality only preserves row and column labels, and does not even
    *attempt* to deal with any supplementary data that could be stored in other
    columns of `anndata.AnnData.obs` or `anndata.AnnData.var`

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

def convert(input_dir: Path) -> Tuple[anndata.AnnData, anndata.AnnData]:
    alevin_dir = input_dir / 'alevin'

    with open(alevin_dir / 'quants_mat_rows.txt') as f:
        cb_names = [line.strip() for line in f]

    with open(alevin_dir / 'quants_mat_cols.txt') as f:
        gene_names = [line.strip() for line in f]

    raw_matrix = scipy.io.mmread(alevin_dir / 'quants_mat.mtx.gz').tocsr()
    raw_labeled = LabeledMatrix(
        matrix=raw_matrix,
        row_labels=cb_names,
        col_labels=gene_names,
    )

    ie_labeled = collapse_intron_cols(raw_labeled)

    density = ie_labeled.matrix.nnz / np.prod(ie_labeled.matrix.shape)
    if density > DENSITY_THRESHOLD:
        ie_labeled.matrix = ie_labeled.matrix.todense()

    return raw_labeled.to_anndata(), ie_labeled.to_anndata()

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('alevin_output_dir', type=Path)
    args = p.parse_args()

    full, summed = convert(args.alevin_output_dir)
    full.write_h5ad('full.h5ad')
    summed.write_h5ad('out.h5ad')
