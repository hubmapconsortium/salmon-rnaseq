#!/usr/bin/env python3
from argparse import ArgumentParser
from functools import lru_cache
from itertools import chain
from pathlib import Path
from subprocess import PIPE, Popen
from typing import Iterable, Tuple

import manhole
from fastq_utils import Read

from common import BARCODE_UMI_FASTQ_PATH, TRANSCRIPT_FASTQ_PATH

SAMTOOLS_COMMAND_TEMPLATE = [
    "samtools",
    "view",
    "{bam_file}",
]


@lru_cache(maxsize=None)
def dummy_barcode_qual(length: int) -> str:
    return "F" * length


def read_bam_file(bam_file: Path) -> Iterable[Tuple[Read, Read]]:
    """
    Convert a BAM file with a call to the `samtools` executable.

    Not using pysam because this doesn't work properly under PyPy,
    and performing barcode adjustment/annotation really benefits
    from that acceleration -- the Docker image for this step can't
    be dynamically selected
    """
    print("Extracting reads from", bam_file)
    command = [piece.format(bam_file=bam_file) for piece in SAMTOOLS_COMMAND_TEMPLATE]
    proc = Popen(command, bufsize=1, stdout=PIPE, universal_newlines=True)
    for line in proc.stdout:
        pieces = line.strip().split("\t")
        read_id = pieces[0]
        seq = pieces[9]
        qual = pieces[10]
        tag_pieces = [p.split(":") for p in pieces[11:]]
        tags = {t[0]: t[2] for t in tag_pieces}
        barcode_umi = tags["XC"] + tags["XM"]

        barcode_umi_read = Read(
            read_id=f"@{read_id}",
            seq=barcode_umi,
            unused="+",
            qual=dummy_barcode_qual(len(barcode_umi)),
        )
        transcript_read = Read(read_id=f"@{read_id}", seq=seq, unused="+", qual=qual)

        yield barcode_umi_read, transcript_read


def main(source_dirs: Iterable[Path], output_dir: Path = Path()):
    buf = output_dir / BARCODE_UMI_FASTQ_PATH
    trf = output_dir / TRANSCRIPT_FASTQ_PATH

    bam_files = chain.from_iterable(d.glob("**/*.bam") for d in source_dirs)

    with open(buf, "w") as buo, open(trf, "w") as tro:
        for bam_file in bam_files:
            for bu, tr in read_bam_file(bam_file):
                print(bu.serialize(), file=buo)
                print(tr.serialize(), file=tro)


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("directory", type=Path)
    args = p.parse_args()

    main(args.directory)
