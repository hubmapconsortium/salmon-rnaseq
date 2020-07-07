#!/usr/bin/env python3
from argparse import ArgumentParser
import json
from os import fspath
from pathlib import Path
from pprint import pprint
from typing import Dict, Iterable, Tuple

import openpyxl
import pandas as pd

FASTQ_EXTENSIONS = [
    'fastq',
    'fastq.gz',
    'fq',
    'fq.gz'
]

FOUND_PAIR_COLOR = '\033[01;32m'
UNPAIRED_COLOR = '\033[01;31m'
NO_COLOR = '\033[00m'

def find_metadata_spreadsheet(directory: Path) -> Path:
    spreadsheets = list(directory.glob('*.xlsx'))
    if len(spreadsheets) != 1:
        addition = ''
        if spreadsheets:
            addition = ':' + ''.join(f'\n\t{s}' for s in spreadsheets)
        message = f'Expected 1 .xlsx file, found {len(spreadsheets)}{addition}'
        raise ValueError(message)
    return spreadsheets[0]

def parse_metadata_spreadsheet(directory: Path, verbose: bool = True) -> Dict[str, str]:
    metadata_xlsx = find_metadata_spreadsheet(directory)
    wb = openpyxl.load_workbook(metadata_xlsx)
    ws = wb['Instrument metadata']
    v = ws.values
    columns = next(v)
    # skip extra description
    next(v)
    df = pd.DataFrame(v, columns=columns)
    row_sel = [any(~row.isnull()) for i, row in df.iterrows()]
    df_usable = df.iloc[row_sel, :]
    # Map R1 FASTQ to dataset ID; we try to find R2 dynamically
    mapping = {
        row.loc['Fastq File - Run 1 - Read1']: row.loc['Data ID']
        for i, row in df_usable.iterrows()
    }
    if verbose:
        print('Metadata mapping:')
        pprint(mapping)
    return mapping

def find_r1_fastq_files(directory: Path) -> Iterable[Path]:
    pattern = '**/*_R1*.{extension}'
    for extension in FASTQ_EXTENSIONS:
        yield from directory.glob(pattern.format(extension=extension))

def find_fastq_files(directory: Path) -> Iterable[Tuple[str, Path, Path]]:
    """
    :param directory:
    :return: Iterable of 3-tuples:
     [0] R1 FASTQ file
     [1] R2 FASTQ file
     [2] dataset label
    """
    fastq_dataset_mapping = parse_metadata_spreadsheet(directory)

    for r1_fastq_file in find_r1_fastq_files(directory):
        dataset_id = fastq_dataset_mapping.get(r1_fastq_file.name)
        if dataset_id is None:
            print(UNPAIRED_COLOR + 'Found FASTQ file not in metadata:' + NO_COLOR)
            print(f'\t{r1_fastq_file}')
            continue

        r2_fastq_filename = r1_fastq_file.name.replace('_R1', '_R2')
        r2_fastq_file = r1_fastq_file.with_name(r2_fastq_filename)

        if r2_fastq_file.is_file():
            print(FOUND_PAIR_COLOR + 'Found pair of FASTQ files:' + NO_COLOR)
            print(f'\t{r1_fastq_file}')
            print(f'\t{r2_fastq_file}')
            yield dataset_id, r1_fastq_file, r2_fastq_file
        else:
            print(UNPAIRED_COLOR + 'Found unpaired FASTQ file:' + NO_COLOR)
            print(f'\t{r1_fastq_file}')

def main(directory: Path):
    sequence_file_bundles = []

    for dataset_id, r1_fastq_path, r2_fastq_path in find_fastq_files(directory):
        sequence_bundle = {
            'fastq_r1': {
                'class': 'File',
                'path': fspath(r1_fastq_path),
            },
            'fastq_r2': {
                'class': 'File',
                'path': fspath(r2_fastq_path),
            },
            'dataset_id': dataset_id,
        }
        sequence_file_bundles.append(sequence_bundle)
    print('Sequence bundles:')
    pprint(sequence_file_bundles)
    output_path = Path('bundles.json')
    print('Saving sequence bundles to', output_path)
    with open(output_path, 'w') as text_file:
        json.dump(sequence_file_bundles, text_file)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('directory', type=Path)
    args = p.parse_args()

    main(args.directory)
