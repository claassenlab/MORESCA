import pytest
import scanpy as sc

from MORESCA.analysis_steps import feature_selection

ADATA = sc.datasets.pbmc3k()
ADATA.layers["counts"] = ADATA.X.copy()


@pytest.mark.parametrize(
    "method", ["seurat", "seurat_v3", "analytical_pearson", "anti_correlation", "triku"]
)
@pytest.mark.parametrize("number_features", [100, 2000, 5000, None])
def test_feature_selection(method, number_features):
    if method == "analytical_pearson" and number_features is None:
        pytest.skip()
    adata = ADATA.copy()
    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.log1p(adata)
    feature_selection(
        adata=adata,
        apply=True,
        method=method,
        number_features=number_features,
        inplace=True,
    )
