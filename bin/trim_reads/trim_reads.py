#!/usr/bin/env python3
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor, wait
from itertools import chain
from os import fspath
from pathlib import Path
from shlex import quote
from shutil import copy
from subprocess import check_call
from typing import Iterable, Sequence, Tuple

from fastq_utils import find_grouped_fastq_files

from common import (
    BARCODE_UMI_FASTQ_PATH,
    TRANSCRIPT_FASTQ_GZ_PATH,
    TRANSCRIPT_FASTQ_PATH,
    Assay,
)

OUTPUT_PATH = Path("trimmed")

TRIM_COMMAND = [
    "/opt/seqtk",
    "trimfq",
    "{input_fastq}",
]


def find_adj_fastq_files(directory: Path) -> Tuple[Path, Path]:
    # not general enough to implement in fastq-utils; very specific
    # to how we create "synthetic" barcode + UMI FASTQ files
    barcode_umi_fastq = directory / BARCODE_UMI_FASTQ_PATH

    transcript_fastq = directory / TRANSCRIPT_FASTQ_PATH
    transcript_fastq_gz = directory / TRANSCRIPT_FASTQ_GZ_PATH

    if transcript_fastq.is_file():
        return barcode_umi_fastq, transcript_fastq
    elif transcript_fastq_gz.is_file():
        return barcode_umi_fastq, transcript_fastq_gz
    else:
        message = (
            f"Couldn't find {TRANSCRIPT_FASTQ_PATH} or {TRANSCRIPT_FASTQ_GZ_PATH} in {directory}"
        )
        raise ValueError(message)


def trim_reads(fastq_r1: Path, fastq_r2: Path, output_subdir: Path):
    print("Copying", quote(fspath(fastq_r1)), "to", quote(fspath(output_subdir)))
    copy(fastq_r1, output_subdir)

    command = [piece.format(input_fastq=fastq_r2) for piece in TRIM_COMMAND]
    fastq_r2_out = output_subdir / fastq_r2.name
    command_str = '"{}"'.format(" ".join(quote(s) for s in command))
    print("Running", command_str, "with output", quote(fspath(fastq_r2_out)))
    with open(fastq_r2_out, "wb") as f:
        check_call(command, stdout=f)


def main(assay: Assay, orig_fastq_dirs: Sequence[Path], adj_fastq_dir: Path, threads: int):
    fastq_pairs: Iterable[Sequence[Path]]
    if assay.barcode_adj_performed:
        if assay.barcode_adj_r1_r2:
            fastq_pairs = find_grouped_fastq_files(adj_fastq_dir, 2)
        else:
            fastq_pairs = [find_adj_fastq_files(adj_fastq_dir)]
    else:
        fastq_pairs = chain.from_iterable(
            find_grouped_fastq_files(fastq_dir, 2) for fastq_dir in orig_fastq_dirs
        )

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = []
        for i, (r1_fastq_file, r2_fastq_file) in enumerate(fastq_pairs, 1):
            subdir = OUTPUT_PATH / str(i)
            subdir.mkdir(exist_ok=True, parents=True)
            future = executor.submit(
                trim_reads,
                r1_fastq_file,
                r2_fastq_file,
                subdir,
            )
            futures.append(future)
        wait(futures)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("adj_fastq_dir", type=Path)
    p.add_argument("orig_fastq_dir", type=Path, nargs="+")
    p.add_argument("-p", "--threads", type=int)
    args = p.parse_args()

    main(args.assay, args.orig_fastq_dir, args.adj_fastq_dir, args.threads)
