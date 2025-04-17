import pytest
import scanpy as sc

from MORESCA.pipeline import clustering

ADATA = sc.datasets.pbmc3k()
ADATA.layers["counts"] = ADATA.X.copy()

sc.pp.filter_cells(ADATA, min_genes=50)
sc.pp.filter_genes(ADATA, min_cells=10)
sc.pp.log1p(ADATA)
sc.pp.pca(ADATA)
sc.pp.neighbors(ADATA)


@pytest.mark.parametrize("method", ["leiden", "phenograph"])
@pytest.mark.parametrize("resolution", [0.5, 1, 2, [0.1, 0.2, 0.3], "auto"])
def test_clustering(method, resolution):
    adata = ADATA.copy()
    clustering(
        adata=adata, resolution=resolution, apply=True, method=method, inplace=True
    )
