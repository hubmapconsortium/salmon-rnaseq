#!/usr/bin/env python3
import csv
import re
from argparse import ArgumentParser
from os import environ, fspath
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

metadata_filename_pattern = re.compile(r"^[0-9A-Fa-f]{32}/probe_set.csv$")

SALMON_COMMAND = [
    "salmon",
    "index",
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
    for file_path in fastq_dir.iterdir():
        if metadata_filename_pattern.match(file_path.name):
            return file_path

def format_index_input_file(index_input_file):
    df = pd.read_csv(index_input_file)
    df = df[["gene_id", "probe_sequence"]]
    formatted_output_file_name = "visium_index_input.csv"
    df.to_csv(formatted_output_file_name, header=False)
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
