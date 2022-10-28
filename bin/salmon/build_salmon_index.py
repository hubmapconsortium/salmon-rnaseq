#!/usr/bin/env python3
import csv
import re
from argparse import ArgumentParser
from os import environ, fspath, walk
from pathlib import Path
from subprocess import check_call
from typing import Iterable, Optional, Sequence, Tuple

import manhole
import pandas as pd

from common import (
    BARCODE_UMI_FASTQ_PATH,
    TRANSCRIPT_FASTQ_GZ_PATH,
    TRANSCRIPT_FASTQ_PATH,
    Assay,
)

probe_set_pattern = "probe_set.csv"

SALMON_COMMAND = [
    "salmon",
    "index",
    "--features",
    "-t",
    "{index_input_file}",
    "-k7",
    "-i",
    "visium_index"
]

def find_index_input_file(fastq_dir):
    """
    Finds and returns the first index input file for a HuBMAP data set.
    """
    for dirpath_str, dirnames, filenames in walk(fastq_dir):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            if filepath.match(probe_set_pattern):
                return filepath

def make_gene_ids_unique(df):
    gene_counts = {}
    for i in df.index:
        gene_id = df.at[i, 'gene_id']
        if gene_id not in gene_counts:
            gene_counts[gene_id] = 0
        else:
            df.at[i, 'gene_id'] = gene_id + '-' + str(gene_counts[gene_id])
        gene_counts[gene_id] += 1
    return df

def format_index_input_file(index_input_file):
    df = pd.read_csv(index_input_file, skiprows=5) #Ignore the comments at the beginning of the file
    df = df[["gene_id", "probe_seq"]]
    df = make_gene_ids_unique(df)
    formatted_output_file_name = "visium_index_input.tsv"
    df.to_csv(formatted_output_file_name, header=False, sep='\t', index=False)
    return formatted_output_file_name

def main(
    assay: Assay,
    fastq_dir: Path,
):

    if assay in {Assay.VISIUM_FFPE}:
        index_input_file = find_index_input_file(fastq_dir)
        formatted_index_input_file = format_index_input_file(index_input_file)
        command = [
            piece.format(
                index_input_file=formatted_index_input_file
            )
            for piece in SALMON_COMMAND
        ]

        print("Running:", " ".join(str(x) for x in command))
        env = environ.copy()
        # Necessary for Singularity; this environment variable isn't
        # set by that container runtime but is required to run Salmon
        env["LD_LIBRARY_PATH"] = "/usr/local/lib"
        check_call(command)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("fastq_dir", type=Path)
    args = p.parse_args()

    main(
        args.assay,
        args.fastq_dir,
    )
