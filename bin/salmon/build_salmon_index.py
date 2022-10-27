#!/usr/bin/env python3
import csv
import re
from argparse import ArgumentParser
from os import environ, fspath
from pathlib import Path
from subprocess import check_call
from typing import Iterable, Optional, Sequence, Tuple

import manhole
from fastq_utils import find_grouped_fastq_files

from common import (
    BARCODE_UMI_FASTQ_PATH,
    TRANSCRIPT_FASTQ_GZ_PATH,
    TRANSCRIPT_FASTQ_PATH,
    Assay,
)

visium_index = ''

SALMON_COMMAND = [
    "salmon",
    "index",
    "-t",
    "{index_input_file}",
    "-i",
    "transcripts_index",
    "-k",
    "31"
]

def find_index_input_file(fastq_dir):
    pass

def format_index_input_file(index_input_file):
    pass

def main(
    assay: Assay,
    fastq_dir: Path,
):

    if assay in {Assay.VISIUM_FFPE}:
        index_input_file = find_index_input_file(fastq_dir)
        formatted_index_input_file = format_index_input_file(index_input_file)
        command = [
            piece.format(
                index_input_file=index_input_file
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
