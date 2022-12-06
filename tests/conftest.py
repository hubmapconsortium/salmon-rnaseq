from pathlib import Path

import pytest


class PathFixture:
    def __init__(self, tmp_path_factory):
        self.tmp_path = tmp_path_factory.mktemp("salmon_rnaseq_fixture")
        self.base_datadir = Path(__file__).parent / "data"

    def get_dir_10x_v2_sn(self):
        return str(self.base_datadir / "10x_v2_sn")


@pytest.fixture(scope="session", autouse=True)
def fixture_setup(tmp_path_factory: pytest.TempPathFactory):
    return PathFixture(tmp_path_factory)
