#!/usr/bin/env python3
"""
Expands a single sciRNA-seq FASTQ file into two: one containing
the original transcript and one with barcode/UMI nucleotide sequences
constructed from the mappings in this repository, for use as input
to Salmon.
"""
import csv
import json
import re
from copy import deepcopy
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, Tuple

from fastq_utils import Read, fastq_reader, smart_open

from common import decompress_fastq

# Relative to path *in container*
BARCODE_DATA_DIR = Path(__file__).parent / "data/sciseq"
FASTQ_INPUT_PATTERN = re.compile(r"(?P<basename>.+)\.(fastq|fq)(.gz)?")

BASE_OUTPUT_DIR = Path("adj_fastq")
METADATA_JSON_PATH = Path("metadata.json")


@lru_cache(maxsize=None)
def get_base_qual_str(length: int) -> str:
    return "F" * length


class BarcodeMapper:
    p5_mapping: Dict[str, str]
    p7_mapping: Dict[str, str]
    rt_mapping: Dict[str, str]
    rt2_mapping: Dict[str, str]
    ligation_mapping: Dict[str, str]

    @classmethod
    def read_barcode_mapping(cls, filename) -> Dict[str, str]:
        mapping = {}
        with open(BARCODE_DATA_DIR / filename) as f:
            for line in f:
                pieces = line.strip().split()
                mapping[pieces[0]] = pieces[1]
        return mapping

    def __init__(self):
        labels = ["p5", "p7", "rt", "rt2", "ligation"]
        for label in labels:
            setattr(self, f"{label}_mapping", self.read_barcode_mapping(f"{label}.txt"))


def convert(mapper: BarcodeMapper, input_fastq: Path, output_dir: Path, basename: str):
    output_dir.mkdir(exist_ok=True, parents=True)
    print("Converting", input_fastq)
    barcode_umi_path = output_dir / f"{basename}_R1.fastq"
    transcript_path = output_dir / f"{basename}_R2.fastq"
    decompress_fastq(input_fastq, transcript_path)

    with smart_open(barcode_umi_path, "wt") as f:
        for transcript_read in fastq_reader(input_fastq):
            id_pieces = transcript_read.read_id.lstrip("@").split("|")
            p7 = mapper.p7_mapping[id_pieces[2]]
            p5 = mapper.p5_mapping[id_pieces[3]]
            rt2_id, lig_id = id_pieces[4].split("_")
            rt2 = mapper.rt2_mapping[rt2_id]
            lig = mapper.ligation_mapping[lig_id]
            umi = id_pieces[5]

            barcode_umi = umi + p7 + p5 + rt2 + lig
            barcode_umi_read = Read(
                read_id=transcript_read.read_id,
                seq=barcode_umi,
                unused=transcript_read.unused,
                qual=get_base_qual_str(len(barcode_umi)),
            )
            print(barcode_umi_read.serialize(), file=f)


def read_barcode_file(barcode_file: Path) -> Dict[str, str]:
    mapping = {}
    with open(barcode_file, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            mapping[row["RT.Barcode"]] = row["Index.Sequence"]
    return mapping


def main(directories: Iterable[Path], output_dir):
    mapper = BarcodeMapper()
    experiment_ids = set()
    extra_barcode_mapping = {}
    for directory in directories:
        mapper_copy = deepcopy(mapper)
        for barcode_file in directory.glob("**/*_barcodes.csv"):
            experiment_barcodes = read_barcode_file(barcode_file)
            extra_barcode_mapping.update(experiment_barcodes)
            mapper_copy.rt2_mapping.update(experiment_barcodes)

        for child in directory.iterdir():
            # no walrus operator; PyPy3 is at 3.6 as of writing this
            m = FASTQ_INPUT_PATTERN.match(child.name)
            if m:
                basename = m.group("basename")
                convert(mapper_copy, child, output_dir, basename)
                experiment_ids.add(basename.split("-")[0])

    # TODO: relax this
    assert len(experiment_ids) == 1
    experiment_id = next(iter(experiment_ids))
    with open(METADATA_JSON_PATH, "w") as f:
        metadata = {
            "experiment_id": experiment_id,
            "extra_barcodes": extra_barcode_mapping,
        }
        json.dump(metadata, f)
