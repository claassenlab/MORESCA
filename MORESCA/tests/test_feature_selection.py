import sys

import pytest
import scanpy as sc

from MORESCA.pipeline import feature_selection

ADATA = sc.datasets.pbmc3k()
ADATA.layers["counts"] = ADATA.X.copy()


@pytest.mark.parametrize(
    "method",
    [
        "seurat",
        "seurat_v3",
        "analytical_pearson",
        "anti_correlation",
        "triku",
        "hotspot",
    ],
)
@pytest.mark.parametrize("number_features", [2000, None])
def test_feature_selection(method, number_features):
    if method == "analytical_pearson" and number_features is None:
        pytest.skip()
    elif method == "anti_correlation" and sys.version_info.minor >= 13:
        # TODO: anticor-features is not installed properly for Python 3.13.0
        # because ray cannot be installed
        pytest.skip(
            "anti_correlation method is not yet compatible with MacOS 13.0 and Python 3.13.0."
        )
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

    if number_features and method not in ["seurat", "anti_correlation"]:
        selected_number_features = adata.var.highly_variable.sum()
        assert selected_number_features == number_features, (
            f"Incorrect number of features after selection, expected {number_features}, found {selected_number_features}."
        )
