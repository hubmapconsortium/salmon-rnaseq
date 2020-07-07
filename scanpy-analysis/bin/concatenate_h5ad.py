#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import List

import anndata

def main(h5ad_files: List[Path]):
    datasets = [anndata.read_h5ad(h5ad_file) for h5ad_file in h5ad_files]
    for dataset in datasets:
        dataset_id = dataset.uns['dataset_id']
        print('Adjusting cell IDs for dataset', dataset_id)
        dataset.obs.loc[:, 'barcode'] = dataset.obs.index
        dataset.obs.loc[:, 'dataset_id'] = dataset_id
        new_index = [f'{dataset_id}_{barcode}' for barcode in dataset.obs.index]
        dataset.obs.index = new_index

    first, *rest = datasets
    concatenated: anndata.AnnData = first.concatenate(*rest, index_unique=None)

    output_path = Path('out.h5ad')
    print('Saving concatenated data to', output_path)
    concatenated.write_h5ad(output_path)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('h5ad_file', type=Path, nargs='+')
    args = p.parse_args()

    main(args.h5ad_file)
