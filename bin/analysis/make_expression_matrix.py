#!/usr/bin/env python3
"""
Takes directory containing a quant-sf file for each sample
These quant-sf files are effectively TSV files with the following fields
Name, Length, Effective Length, TPM, and NumReads
Creates a hdf file containing matrices for TPM and NumReads
"""
from argparse import ArgumentParser
from os import fspath
from pathlib import Path
from typing import List, Tuple

import manhole
import pandas as pd


def get_sample_id(quant_file_name: str) -> str:
    """
    Get the sample id from the file name
    """
    return quant_file_name.split("/")[-1][:-9]


def initialize_matrix_dfs(
    source_df: pd.DataFrame, sample_ids: List
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    transcript_names = source_df["Name"]

    tpm_df = pd.DataFrame(
        transcript_names,
        columns=sample_ids,
        index=transcript_names,
    )
    num_reads_df = pd.DataFrame(
        transcript_names,
        columns=sample_ids,
        index=transcript_names,
    )

    return tpm_df, num_reads_df


def main(quant_dir: Path):
    initialized = False

    num_reads_df = None
    tpm_df = None

    quant_files = [fspath(quant_file) for quant_file in list(quant_dir.glob("**/*quant.sf"))]

    sample_ids = [get_sample_id(fspath(quant_file)) for quant_file in list(quant_files)]

    # For each sample
    for quant_file in quant_files:
        sample_id = get_sample_id(fspath(quant_file))  # Get the sample_id
        source_df = pd.read_csv(fspath(quant_file), delimiter="\t")  # And open the source file

        if not initialized:  # If we haven't initialized our matrices yet
            tpm_df, num_reads_df = initialize_matrix_dfs(source_df, sample_ids)  # Do it now
            initialized = True

        tpm_df[sample_id] = source_df["TPM"].values
        num_reads_df[sample_id] = source_df["NumReads"].values

    tpm_df = tpm_df.transpose()
    num_reads_df = num_reads_df.transpose()

    # Write out to hdf5 file
    tpm_df.to_hdf("expression_matrices.h5", key="tpm", mode="w")
    num_reads_df.to_hdf("expression_matrices.h5", key="num_reads")


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    # Make this function take a directory and get the quant_sf files from it instead
    p.add_argument("quant_dir", type=Path)
    args = p.parse_args()

    main(args.quant_dir)
