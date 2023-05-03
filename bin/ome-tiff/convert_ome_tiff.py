#!/usr/bin/env python3
import json
import re
import shlex
from argparse import ArgumentParser
from collections import defaultdict
from os import walk
from pathlib import Path
from subprocess import run
from typing import Iterable, Tuple
from aicsimageio import AICSImage

ome_tiff_pattern = re.compile(r"(?P<basename>.*)\.tif(f?)$")

def find_ome_tiffs(input_dir: Path) -> Iterable[Tuple[Path, Path]]:
    """
    Yields 2-tuples:
     [0] full Path to source file
     [1] output file Path (source file relative to input_dir)
    """
    for dirpath_str, _, filenames in walk(input_dir):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            if ome_tiff_pattern.match(filename):
                src_filepath = dirpath / filename
                dest_filepath = dirpath / filename.replace(".tiff", ".ome.tiff")
                yield src_filepath, dest_filepath.relative_to(input_dir)


def fix_ome_tiff(source: Path, dest: Path):
    img = AICSImage(source)
    img.write(dest)


def main(input_dir: Path, output_path_prefix):
    for source, dest_relative in find_ome_tiffs(input_dir):
        dest_absolute = output_path_prefix / dest_relative
        fix_ome_tiff(source, dest_absolute)

if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("input_dir", type=Path)
    p.add_argument("--output-path-prefix", type=Path, default=Path())
    args = p.parse_args()

    main(args.input_dir, args.output_path_prefix)