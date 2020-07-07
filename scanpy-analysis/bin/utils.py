from enum import Enum
import bz2
import gzip
import lzma
from os import PathLike
from pathlib import Path

class FileType(Enum):
    def __new__(cls, filetype, open_function):
        obj = object.__new__(cls)
        obj._value_ = filetype
        obj.open_function = open_function
        return obj

    GZ = ("gz", gzip.open)
    BZ2 = ("bz2", bz2.open)
    XZ = ("xz", lzma.open)
    TEXT = ("txt", open)

def get_file_type_by_extension(file_path: Path) -> FileType:
    suffix = file_path.suffix.lstrip(".")
    try:
        return FileType(suffix)
    except ValueError:
        # No special suffix, assume text
        return FileType.TEXT

def smart_open(file_path: PathLike, mode="rt", *args, **kwargs):
    file_type = get_file_type_by_extension(Path(file_path))
    return file_type.open_function(file_path, mode, *args, **kwargs)
