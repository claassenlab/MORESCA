import os
from pathlib import Path
from typing import List, Literal, Tuple, Union

import numpy as np
import pytest
import scanpy as sc

from MORESCA.pipeline import (
    clustering,
    feature_selection,
    load_data,
    neighborhood_graph,
    pca,
)

ADATA = sc.datasets.pbmc3k()
ADATA.layers["counts"] = ADATA.X.copy()
DATA_PATH_PREFIX = "data/adata"


@pytest.mark.parametrize(
    "ext, as_path", [("h5ad", True), ("h5ad", False), ("loom", False)]
)
def test_load_data(ext: str, as_path: bool):
    DATA_PATH = f"{DATA_PATH_PREFIX}.{ext}"
    if ext == "loom":
        ADATA.write_loom(DATA_PATH)
    else:
        ADATA.write_h5ad(DATA_PATH)
    if as_path:
        DATA_PATH = Path(DATA_PATH)
    adata = load_data(DATA_PATH)
    assert isinstance(adata, sc.AnnData), "Return value should be an AnnData object."
    os.remove(DATA_PATH)


@pytest.mark.parametrize(
    "apply, inplace", [(True, True), (True, False), (False, False)]
)
@pytest.mark.parametrize(
    "resolution, run_pca",
    [
        (0.5, True),
        (1, True),
        (0.5, False),
        (1, False),
        ([0.5, 1], True),
        ((0.5, 1), True),
        ("auto", True),
    ],
)
def test_clustering(
    apply: bool,
    resolution: Union[
        float, int, List[Union[float, int]], Tuple[Union[float, int]], Literal["auto"]
    ],
    inplace: bool,
    run_pca: bool,
):
    adata = ADATA.copy()
    feature_selection(
        adata=adata, apply=True, method="seurat_v3", number_features=1000, inplace=True
    )
    if run_pca:
        pca(adata=adata, apply=True, inplace=True)
    neighborhood_graph(adata=adata, apply=True, inplace=True)
    return_val = clustering(
        adata=adata, apply=apply, resolution=resolution, inplace=inplace
    )

    # Check correct return type
    if inplace or not apply:
        assert return_val is None, "Return value should be None if inplace=True"
        if not apply:
            assert not adata.obs.columns.str.startswith("leiden").any(), (
                "Clustering should not have been applied."
            )
            return
        adata_ = adata
    else:
        assert isinstance(return_val, sc.AnnData), (
            "Return value should be an AnnData object"
        )
        adata_ = return_val

    # Check if the clustering is stored in the obs DataFrame
    if isinstance(resolution, (int, float)):
        assert f"leiden_r{resolution}" in adata_.obs.columns
    elif isinstance(resolution, (list, tuple)) or resolution == "auto":
        resolutions = (
            resolution
            if resolution != "auto"
            else adata_.uns["MORESCA"]["clustering"]["resolutions"]
        )
        for res in resolutions:
            assert f"leiden_r{res}" in adata_.obs.columns, (
                f"`leiden_r{res}` should be in obs DataFrame."
            )
        assert "leiden" in adata_.obs.columns, "`leiden` should be in obs DataFrame."
        assert "silhouette_scores" in adata_.uns["MORESCA"]["clustering"].keys(), (
            "Scores should be stored."
        )
        assert len(adata_.uns["MORESCA"]["clustering"]["silhouette_scores"]) == len(
            resolutions
        ), "Scores should be stored for each resolution."
        assert np.all(
            adata_.obs["leiden"]
            == adata_.obs[
                f"leiden_r{resolutions[np.argmax(adata_.uns['MORESCA']['clustering']['silhouette_scores'])]}"
            ]
        ), "leiden should be the leiden_r with the highest score."


def test_clustering_exceptions():
    pass
