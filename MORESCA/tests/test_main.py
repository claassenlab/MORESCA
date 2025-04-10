import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import scanpy as sc

from MORESCA.main import run_analysis

ADATA = sc.datasets.pbmc3k()
ADATA.obs["batch"] = 1350 * ["a"] + 1350 * ["b"]
DATA_PATH = Path("data/data_raw.h5ad")
ADATA.write_h5ad(DATA_PATH)
CONFIG_PATH = Path("test-config.gin")


def test_run_analysis():
    run_analysis(DATA_PATH, config_path=CONFIG_PATH)


@pytest.mark.parametrize(
    "data_path, output_path",
    [
        ([DATA_PATH], ["res/1", "res/2"]),
        ([DATA_PATH, DATA_PATH, DATA_PATH], ["res/1", "res/2"]),
    ],
)
def test_run_analysis_exception(data_path: list, output_path: list):
    with pytest.raises(
        ValueError, match="Incompatible values for `result_path` and `data_path`"
    ):
        run_analysis(
            data_path=data_path, config_path=CONFIG_PATH, result_path=output_path
        )


@pytest.mark.parametrize(
    "data_path, output_path",
    [
        ([DATA_PATH], ["res"]),
        ([DATA_PATH, DATA_PATH], ["res"]),
        ([DATA_PATH, DATA_PATH], ["res/1", "res/2"]),
    ],
)
def test_main(data_path: list, output_path: list):
    cmd = (
        [sys.executable, "../main.py", "-d"]
        + data_path
        + ["-o"]
        + output_path
        + ["-p", CONFIG_PATH]
    )
    subprocess.run(cmd, check=True)
    sample_str = ""
    for i, path in enumerate(output_path):
        if len(data_path) > 1 and len(output_path) == 1:
            sample_str = f"_s{i + 1:02d}"
        assert os.path.isfile(os.path.join(path, f"data_processed{sample_str}.h5ad")), (
            "Output file not found."
        )
    shutil.rmtree("res")


def test_moresca():
    cmd = ["moresca", "-d", DATA_PATH, "-p", CONFIG_PATH, "-v"]
    subprocess.run(cmd, check=True)
