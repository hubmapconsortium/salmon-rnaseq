#!/usr/bin/env python3
import csv
import re
from argparse import ArgumentParser
from os import environ, fspath
from pathlib import Path
from subprocess import check_call
from typing import Iterable, Optional, Sequence, Tuple

from fastq_utils import find_grouped_fastq_files

from common import (
    BARCODE_UMI_FASTQ_PATH,
    TRANSCRIPT_FASTQ_GZ_PATH,
    TRANSCRIPT_FASTQ_PATH,
    Assay,
)

SALMON_COMMAND = [
    "salmon",
    "alevin",
    "--index",
    "/opt/gencode.v35.intron-exon.sidx",
    "--libType",
    "A",
    "--output",
    "salmon_out",
    "--dumpMtx",
    "{salmon_option}",
    "--tgMap",
    "/opt/gencode.v35.annotation.expanded.tx2gene.tsv",
    "-p",
    "{threads}",
]

cell_count_filename = "extras/expected_cell_count.txt"
metadata_filename_pattern = re.compile(r"^[0-9A-Fa-f]{32}-metadata.tsv$")
metadata_cell_count_field = "expected_cell_count"


def find_metadata_file(directory: Path) -> Optional[Path]:
    """
    Finds and returns the first metadata file for a HuBMAP data set.
    Does not check whether the dataset ID (32 hex characters) matches
    the directory name, nor whether there might be multiple metadata files.
    """
    for file_path in directory.iterdir():
        if metadata_filename_pattern.match(file_path.name):
            return file_path


def read_expected_cell_count(directory: Path) -> Optional[int]:
    cell_count_from_file = None
    cell_count_metadata = None

    cell_count_file = directory / cell_count_filename
    if cell_count_file.is_file():
        with open(cell_count_file) as f:
            cell_count_from_file = int(f.read().strip())
            print(f"Read expected cell count from {cell_count_file}: {cell_count_from_file}")

    maybe_metadata_file = find_metadata_file(directory)
    if maybe_metadata_file.is_file():
        with open(maybe_metadata_file, newline="") as f:
            r = csv.DictReader(f, delimiter="\t")
            metadata = next(r)
            if (
                metadata_cell_count_field in metadata
                and metadata[metadata_cell_count_field].isdigit()
            ):
                cell_count_metadata = int(metadata[metadata_cell_count_field])
                print(
                    f"Read expected cell count from {maybe_metadata_file}: {cell_count_metadata}"
                )

    present_cell_counts = sum(x is not None for x in [cell_count_from_file, cell_count_metadata])
    if present_cell_counts == 0:
        return None
    elif present_cell_counts == 1:
        return cell_count_from_file or cell_count_metadata
    else:
        if cell_count_from_file == cell_count_metadata:
            return cell_count_from_file
        else:
            message = (
                f"Found mismatched cell counts: {cell_count_from_file} in {cell_count_file}, "
                f"and {cell_count_metadata} in {maybe_metadata_file}"
            )
            raise ValueError(message)


def read_expected_cell_counts(directories: Sequence[Path]) -> Optional[int]:
    cell_counts = []
    for directory in directories:
        cell_count = read_expected_cell_count(directory)
        if cell_count is not None:
            cell_counts.append(cell_count)

    dirs_with_cell_counts = len(cell_counts)
    if dirs_with_cell_counts == 0:
        return None
    elif dirs_with_cell_counts == len(directories):
        total_expected = sum(cell_counts)
        print("Total expected cells:", total_expected)
        return total_expected
    else:
        message = (
            f"Found expected cell counts in {dirs_with_cell_counts} of "
            f"{len(directories)} directories, need 0 or {len(directories)} "
            f"input directories with cell counts (can't mix auto-detection "
            f"and guided cell barcode counting)"
        )
        raise ValueError(message)


def find_adj_fastq_files(directory: Path) -> Iterable[Tuple[Path, Path]]:
    # not general enough to implement in fastq-utils; very specific
    # to how we create "synthetic" barcode + UMI FASTQ files
    for subdir in directory.iterdir():
        barcode_umi_fastq = subdir / BARCODE_UMI_FASTQ_PATH

        transcript_fastq = subdir / TRANSCRIPT_FASTQ_PATH
        transcript_fastq_gz = subdir / TRANSCRIPT_FASTQ_GZ_PATH

        if transcript_fastq.is_file():
            yield barcode_umi_fastq, transcript_fastq
        elif transcript_fastq_gz.is_file():
            yield barcode_umi_fastq, transcript_fastq_gz


def find_slideseq_barcode_file(base_dir: Path) -> Path:
    pattern = "**/*matched_bead_barcodes.txt"
    barcode_files = list(base_dir.glob(pattern))
    if len(barcode_files) != 1:
        message_pieces = [
            f"Need exactly 1 file matching {pattern} "
            f"under {base_dir}, found {len(barcode_files)}:"
        ]
        message_pieces.extend(f"\t{bf}" for bf in barcode_files)
        raise ValueError("\n".join(message_pieces))
    return barcode_files[0]


def main(
    assay: Assay,
    orig_fastq_dirs: Sequence[Path],
    trimmed_fastq_dir: Path,
    expected_cell_count: Optional[int],
    threads: Optional[int],
):
    threads = threads or 1
    command = [
        piece.format(
            salmon_option=assay.salmon_option,
            threads=threads,
        )
        for piece in SALMON_COMMAND
    ]

    fastq_pairs: Iterable[Sequence[Path]]
    if assay.barcode_adj_performed:
        if assay.barcode_adj_r1_r2:
            fastq_pairs = list(find_grouped_fastq_files(trimmed_fastq_dir, 2))
        else:
            fastq_pairs = list(find_adj_fastq_files(trimmed_fastq_dir))
    else:
        fastq_pairs = list(find_grouped_fastq_files(trimmed_fastq_dir, 2))

    if not fastq_pairs:
        raise ValueError("No FASTQ files found")

    if assay.keep_all_barcodes:
        command.extend(["--keepCBFraction", "1"])
    # hack
    if assay == Assay.SLIDESEQ:
        # Don't support multiple input directories for Slide-seq; this will
        # likely cause significantly incorrect results due to barcode overlap
        # between multiple input data sets
        if len(orig_fastq_dirs) != 1:
            raise ValueError("Need exactly 1 input directory for Slide-seq")
        barcode_file = find_slideseq_barcode_file(orig_fastq_dirs[0])
        command.extend(["--whitelist", fspath(barcode_file)])

    maybe_cell_count = read_expected_cell_count(orig_fastq_dirs) or expected_cell_count
    if maybe_cell_count is not None:
        command.extend(["--forceCells", str(maybe_cell_count)])

    for r1_fastq_file, r2_fastq_file in fastq_pairs:
        fastq_extension = [
            "-1",
            r1_fastq_file,
            "-2",
            r2_fastq_file,
        ]
        command.extend(fastq_extension)

    print("Running:", " ".join(str(x) for x in command))
    env = environ.copy()
    # Necessary for Singularity; this environment variable isn't
    # set by that container runtime but is required to run Salmon
    env["LD_LIBRARY_PATH"] = "/usr/local/lib"
    check_call(command)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("trimmed_fastq_dir", type=Path)
    p.add_argument("orig_fastq_dir", type=Path, nargs="+")
    p.add_argument("--expected-cell-count", type=int)
    p.add_argument("-p", "--threads", type=int)
    args = p.parse_args()

    main(
        args.assay,
        args.orig_fastq_dir,
        args.trimmed_fastq_dir,
        args.expected_cell_count,
        args.threads,
    )
