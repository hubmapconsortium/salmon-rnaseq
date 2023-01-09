#!/usr/bin/env python3
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor, wait
from os import fspath
from pathlib import Path
from subprocess import check_call
from typing import Tuple

import manhole
from fastq_utils import collect_fastq_files_by_directory

FASTQC_COMMAND_TEMPLATE = [
    "fastqc",
    "--outdir",
    "{out_dir}",
]


def single_file_fastqc(fastq_file_and_subdir: Tuple[Path, Path]):
    """
    Run FastQC on a single fastq file

    Takes an absolute path to the input file and a relative path to
    the output subdirectory
    """
    command = [piece.format(out_dir=fastq_file_and_subdir[1]) for piece in FASTQC_COMMAND_TEMPLATE]
    command.append(fspath(fastq_file_and_subdir[0]))
    print("Running:", " ".join(command))
    check_call(command)


def main(directory: Path, threads: int):
    """
    Crawl directory, create appropriate output subdirectories based on input directory structure
    Append output files to a list to pass to Pool.imap_unordered
    """
    fastq_files_by_directory = collect_fastq_files_by_directory(directory)
    print("Found", len(fastq_files_by_directory), "directories containing FASTQ files")

    fastqc_out_dir = Path("fastqc_output")

    fastq_files_and_subdirs = []

    for directory, files in fastq_files_by_directory.items():
        subdir = fastqc_out_dir / directory
        subdir.mkdir(exist_ok=True, parents=True)
        for file in files:
            fastq_files_and_subdirs.append((file, subdir))

    with ProcessPoolExecutor(max_workers=threads) as pool:
        futures = {
            pool.submit(single_file_fastqc, fastq_file_and_subdir): fastq_file_and_subdir
            for fastq_file_and_subdir in fastq_files_and_subdirs
        }
        wait(futures)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("directory", type=Path)
    p.add_argument("threads", type=int)
    args = p.parse_args()

    main(args.directory, args.threads)
