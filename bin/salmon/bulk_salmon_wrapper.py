#!/usr/bin/env python3
from argparse import ArgumentParser
from os import fspath
from pathlib import Path
from subprocess import check_call
from typing import Optional

import manhole
from fastq_utils import find_grouped_fastq_files, get_sample_id_from_r1

human_index = "/opt/gencode.v35.intron-exon.sidx"
human_transcript_map = "/opt/gencode.v35.annotation.expanded.tx2gene.tsv"
mouse_index = "/opt/gencode.vM28.intron-exon.sidx"
mouse_transcript_map = "/opt/gencode.vM28.annotation.expanded.tx2gene.tsv"

SALMON_COMMAND = [
    "salmon",
    "quant",
    "--index",
    "{index}",
    "--libType",
    "A",
    "-p",
    "{threads}",
    "--output",
    "out",
]

TAR_AND_ZIP_COMMAND = [
    "tar",
    "-czvf",
    "out/aux_files.tar.gz",
    "out/aux_info",
]


def rename_file(old_file_name: str, new_file_name: str):
    command = ["mv"]
    command.append(old_file_name)
    command.append(new_file_name)
    check_call(command)


def main(threads: int, directory: Path, organism: Optional[str] = "human"):
    index = human_index if organism == "human" else mouse_index
    for r1_fastq_file, r2_fastq_file in find_grouped_fastq_files(directory, 2):
        command = [piece.format(threads=threads, index=index) for piece in SALMON_COMMAND]

        fastq_extension = [
            "-1",
            fspath(r1_fastq_file),
            "-2",
            fspath(r2_fastq_file),
        ]

        command.extend(fastq_extension)
        print("Running:", command)
        check_call(command)

        check_call(TAR_AND_ZIP_COMMAND)
        # tar and zip auxiliary files

        sample_id = get_sample_id_from_r1(r1_fastq_file)

        # Tag output files with sample_id
        rename_file("out/quant.sf", f"out/{sample_id}-quant.sf")
        rename_file("out/cmd_info.json", f"out/{sample_id}-cmd_info.json")
        rename_file("out/aux_files.tar.gz", f"out/{sample_id}-aux_files.tar.gz")


if __name__ == "__main__":
    manhole.install(activate_on="USR1")

    p = ArgumentParser()
    p.add_argument("-p", "--threads", type=int)
    p.add_argument("directory", type=Path)
    p.add_argument("--organism", type=str, nargs="?", default="human")
    args = p.parse_args()

    main(args.threads, args.directory, args.organism)
