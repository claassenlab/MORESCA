import pytest
import scanpy as sc

from MORESCA.pipeline import normalization

ADATA = sc.datasets.pbmc3k()
ADATA.layers["counts"] = ADATA.X.copy()


@pytest.mark.parametrize("remove_mt", [True, False])
@pytest.mark.parametrize("remove_rb", [True, False])
@pytest.mark.parametrize("remove_hb", [True, False])
@pytest.mark.parametrize(
    "method", ["log1pCP10k", "log1pPF", "PFlog1pPF", "analytical_pearson"]
)
def test_normalization(method, remove_mt, remove_rb, remove_hb):
    adata = ADATA.copy()
    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_genes(adata, min_cells=10)
    normalization(
        adata=adata,
        apply=True,
        method=method,
        inplace=True,
        remove_mt=remove_mt,
        remove_rb=remove_rb,
        remove_hb=remove_hb,
    )
