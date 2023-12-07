import os.path

import gin
import scanpy as sc
import pytest

from MORESCA.analysis_steps import (
    batch_effect_correction,
    feature_selection,
    normalization,
    pca,
    quality_control,
    scaling,
    umap,
    plotting,
    neighborhood_graph,
)

# PATH_H5AD = Path("../data/data_raw.h5ad")
# STR_H5AD = "../data/data_raw.h5ad"
ADATA_original = sc.datasets.pbmc3k()
# sc.read(PATH_H5AD)
ADATA_original.layers["counts"] = ADATA_original.X.copy()

ADATA = None

gin.parse_config_file("test-config.gin")


# Todo: Adapt the CI tests to automatically download the data.
# @pytest.mark.parametrize("test_data_path", [(PATH_H5AD), (STR_H5AD)])
# def test_load_data(test_data_path):
#    assert load_data(test_data_path)


# SOURCE: https://stackoverflow.com/a/22638709
@pytest.fixture(autouse=True)
def run_before_tests():
    global ADATA
    ADATA = ADATA_original.copy()


def test_quality_control():
    quality_control(adata=ADATA)


def test_normalization():
    normalization(adata=ADATA)


def test_feature_selection():
    feature_selection(adata=ADATA)


def test_scaling():
    scaling(adata=ADATA)


def test_pca():
    feature_selection(adata=ADATA)
    pca(adata=ADATA)


def test_umap_with_pca():
    feature_selection(adata=ADATA)
    neighborhood_graph(adata=ADATA)
    pca(adata=ADATA)
    umap(adata=ADATA, pca_before_umap=True)
    assert "neighbors" in ADATA.uns


def test_umap_without_pca():
    umap(adata=ADATA, pca_before_umap=False)
    assert "neighbors_without_pca" in ADATA.uns


def test_plotting_umap():
    umap(adata=ADATA, pca_before_umap=False)
    plotting(adata=ADATA, umap=True, path="figures/")
    assert os.path.exists("figures/umap.png")

    # Clean up after test
    os.remove("figures/umap.png")


def test_plotting_umap_with_pca():
    feature_selection(adata=ADATA)
    pca(adata=ADATA)
    umap(adata=ADATA, pca_before_umap=True)
    plotting(adata=ADATA, umap=True, path="figures/")
    assert os.path.exists("figures/umap.png")

    # Clean up after test
    os.remove("figures/umap.png")


def test_batch_effect_correction():
    ADATA_BATCH = ADATA
    ADATA_BATCH.obs["batch"] = 1350 * ["a"] + 1350 * ["b"]

    pca(ADATA_BATCH, use_highly_variable=False)
    batch_effect_correction(adata=ADATA_BATCH)
