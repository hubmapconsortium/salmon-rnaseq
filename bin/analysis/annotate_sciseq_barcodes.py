#!/usr/bin/env python3
import json
from argparse import ArgumentParser
from pathlib import Path
from typing import Dict, Optional

import anndata
import pandas as pd

BARCODE_DATA_DIR = here = Path(__file__).parent / "data/sciseq"

BARCODE_LENGTH = 10
BARCODE_STARTS = [0, 10, 20]
BARCODE_SEGMENTS = [slice(start, start + BARCODE_LENGTH) for start in BARCODE_STARTS]
UMI_SEGMENT = slice(0, 8)


class CellIdMapper:
    p5_mapping: Dict[str, str]
    p7_mapping: Dict[str, str]
    rt_mapping: Dict[str, str]
    rt2_mapping: Dict[str, str]

    @classmethod
    def read_barcode_mapping(cls, filename) -> Dict[str, str]:
        mapping = {}
        with open(BARCODE_DATA_DIR / filename) as f:
            for line in f:
                pieces = line.strip().split()
                mapping[pieces[1]] = pieces[0]
        return mapping

    def __init__(self, extra_barcodes: Optional[Dict[str, str]]):
        labels = ["p5", "p7", "rt", "rt2"]
        for label in labels:
            setattr(self, f"{label}_mapping", self.read_barcode_mapping(f"{label}.txt"))

        if extra_barcodes is not None:
            self.rt2_mapping.update(extra_barcodes)


def annotate(mapper: CellIdMapper, data: anndata.AnnData, experiment_id: str) -> anndata.AnnData:
    cell_ids = []

    for i in data.obs.index:
        barcode_pieces = [i[seg] for seg in BARCODE_SEGMENTS]
        p7 = mapper.p7_mapping[barcode_pieces[0]]
        p5 = mapper.p5_mapping[barcode_pieces[1]]
        rt2 = mapper.rt2_mapping[barcode_pieces[2]]

        cell_id = f"{p7}_{p5}_{rt2}_{experiment_id}"
        cell_ids.append(cell_id)

    data.obs.loc[:, "cell_id"] = pd.Series(cell_ids, index=data.obs.index)

    return data


def main(h5ad_file: Path, metadata_json_file: Path):
    data = anndata.read_h5ad(h5ad_file)
    with open(metadata_json_file) as f:
        metadata = json.load(f)

    experiment_id = metadata["experiment_id"]

    extra_barcodes = None
    if "extra_barcodes" in metadata:
        extra_barcodes = {value: key for key, value in metadata["extra_barcodes"].items()}

    mapper = CellIdMapper(extra_barcodes)
    return annotate(mapper, data, experiment_id)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("h5ad_file", type=Path)
    p.add_argument("metadata_json_file", type=Path)
    args = p.parse_args()

    d = main(args.h5ad_file, args.metadata_json_file)
    d.write_h5ad("expr.h5ad")
