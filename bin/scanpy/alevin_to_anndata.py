#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

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

def convert(input_dir: Path) -> anndata.AnnData:
    alevin_dir = input_dir / 'alevin'

    with open(alevin_dir / 'quants_mat_rows.txt') as f:
        cb_names = [line.strip() for line in f]
    obs_df = pd.DataFrame(index=cb_names)

    with open(alevin_dir / 'quants_mat_cols.txt') as f:
        gene_names = [line.strip() for line in f]
    var_df = pd.DataFrame(index=gene_names)

    matrix = sparse_matrix = scipy.io.mmread(alevin_dir / 'quants_mat.mtx.gz').tocsr()

    density = sparse_matrix.nnz / np.prod(sparse_matrix.shape)
    if density > DENSITY_THRESHOLD:
        matrix = sparse_matrix.todense()

    d = anndata.AnnData(X=matrix, obs=obs_df, var=var_df)
    return d

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('alevin_output_dir', type=Path)
    args = p.parse_args()

    alv = convert(args.alevin_output_dir)
    alv.write_h5ad('out.h5ad')
