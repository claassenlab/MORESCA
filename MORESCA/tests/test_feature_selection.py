import os
from pathlib import Path
from typing import List, Literal, Tuple, Union

import numpy as np
import pytest
import scanpy as sc

from MORESCA.analysis_steps import (
    clustering,
    feature_selection,
    load_data,
    neighborhood_graph,
    pca,
)

ADATA = sc.datasets.pbmc3k()
ADATA.layers["counts"] = ADATA.X.copy()


@pytest.mark.parametrize(
    "method, number_features",
    [
        ("seurat", 2000),
        ("seurat_v3", 2000),
        ("analytical_pearson", 2000),
        ("anti_correlation", None),
        ("triku", None)
    ],
)
def test_feature_selection(method, number_features):
    adata = ADATA.copy()
    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_genes(adata, min_cells=10)
    sc.pp.log1p(adata)
    feature_selection(
        adata=adata, apply=True, method=method, number_features=number_features, inplace=True
    )
