#!/usr/bin/env python3
"""
Adapted from https://github.com/mruffalo/vpolo/blob/master/vpolo/alevin/parser.py
which is licensed under the GPL v3 (so this file is as well)
"""
import gzip
import sys
from argparse import ArgumentParser
from pathlib import Path
from struct import Struct

import anndata
import numpy as np
import pandas as pd
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

def convert(input_dir: Path, dataset_id: str) -> anndata.AnnData:
    """
    Read the quants sparse binary output of Alevin and converts to an `anndata` object
    """
    data_type = "f"

    alevin_dir = input_dir / 'alevin'

    gene_names = pd.read_csv(alevin_dir / 'quants_mat_cols.txt', header=None)[0].values
    cb_names = pd.read_csv(alevin_dir / 'quants_mat_rows.txt', header=None)[0].values
    num_genes = len(gene_names)
    num_entries = int(np.ceil(num_genes/8))

    obs_df = pd.DataFrame(index=cb_names)
    var_df = pd.DataFrame(index=gene_names)

    with gzip.open(alevin_dir / 'quants_mat.gz') as f:
        line_count = 0
        tot_umi_count = 0

        entries = []
        col_indices = []
        row_indices = []

        header_struct = Struct("B" * num_entries)
        cell_index = 0
        while True:
            gene_index = 0
            line_count += 1
            if not (line_count % 100):
                print("\rDone reading", line_count, "cells", end="")
                sys.stdout.flush()
            try:
                num_exp_genes = 0
                exp_counts = header_struct.unpack_from(f.read(header_struct.size))
                for exp_count in exp_counts:
                    num_exp_genes += bin(exp_count).count("1")

                data_struct = Struct(data_type * num_exp_genes)
                sparse_cell_counts_vec = list(data_struct.unpack_from(f.read(data_struct.size)))[::-1]
                cell_umi_counts = sum(sparse_cell_counts_vec)

            except Exception:
                print("\nRead total", line_count - 1, " cells")
                print("Found total", tot_umi_count, "reads")
                break

            if cell_umi_counts > 0.0:
                tot_umi_count += cell_umi_counts

                cell_counts_vec = []
                for exp_count in exp_counts:
                    for bit in format(exp_count, '08b'):
                        if len(cell_counts_vec) >= num_genes:
                            break

                        if bit == '0':
                            cell_counts_vec.append(0.0)
                        else:
                            abund = sparse_cell_counts_vec.pop()
                            entries.append(abund)
                            col_indices.append(gene_index)
                            row_indices.append(cell_index)
                            cell_counts_vec.append(abund)

                        gene_index += 1

                if len(sparse_cell_counts_vec) > 0:
                    print("Failure in consumption of data")
                    print("left with {} entry(ies)".format(len(sparse_cell_counts_vec)))
            else:
                raise ValueError("Found a CB with no read count, something is wrong")

            cell_index += 1

    matrix = sparse_matrix = scipy.sparse.coo_matrix(
        (entries, (row_indices, col_indices)),
        shape=(len(cb_names), num_genes),
    ).tocsr()
    density = sparse_matrix.nnz / np.prod(sparse_matrix.shape)
    if density > DENSITY_THRESHOLD:
        matrix = sparse_matrix.todense()

    d = anndata.AnnData(X=matrix, obs=obs_df, var=var_df)
    d.uns['dataset_id'] = dataset_id
    return d

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('alevin_output_dir', type=Path)
    p.add_argument('dataset_id')
    args = p.parse_args()

    alv = convert(args.alevin_output_dir, args.dataset_id)
    alv.write_h5ad('out.h5ad')
