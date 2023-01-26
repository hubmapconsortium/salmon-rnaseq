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

def trim_reads_visium_ffpe(fastq_file_path):
    new_sections = []
    text = open(fastq_file_path).read()
    sections = text.split('+')
    for section in sections:
        lines = section.split('\n')
        sequence = lines[2]
        sequence = sequence[0:50]
        lines[2] = sequence
        section = '\n'.join(lines)
        new_sections.append(section)
    new_text = '+'.join(new_sections)
    fastq_file_stem = fastq_file_path.stem
    new_fastq_file_path = Path(f'trimmed_{fastq_file_stem}.fastq')
    new_fastq_file = open(new_fastq_file_path)
    new_fastq_file.write(new_text)
    return

def main(assay, orig_fastq_dirs: Sequence[Path], adj_fastq_dir: Path, threads: int):
    fastq_pairs = chain.from_iterable(
        find_grouped_fastq_files(fastq_dir, 2) for fastq_dir in orig_fastq_dirs
    )

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = []
        for i, (r1_fastq_file, r2_fastq_file) in enumerate(fastq_pairs, 1):
            subdir = OUTPUT_PATH / str(i)
            subdir.mkdir(exist_ok=True, parents=True)
            future = executor.submit(
                trim_reads_visium_ffpe,
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
