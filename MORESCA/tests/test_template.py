import scanpy as sc
from pathlib import Path
from MORESCA.template import run_analysis

ADATA = sc.datasets.pbmc3k()
ADATA.obs["batch"] = 1350 * ["a"] + 1350 * ["b"]
DATA_PATH = Path("data/data_raw.h5ad")
ADATA.write_h5ad(DATA_PATH)


def test_run_analysis():
    run_analysis(DATA_PATH, config_path=Path("test-config.gin"), figures=False, verbose=False)
