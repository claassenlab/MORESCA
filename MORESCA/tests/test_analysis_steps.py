from pathlib import Path

import gin
import pytest
import scanpy as sc

from MORESCA.analysis_steps import (
    batch_effect_correction,
    feature_selection,
    load_data,
    normalization,
    pca,
    quality_control,
    scaling,
)

PATH_H5AD = Path("../data/data_raw.h5ad")
STR_H5AD = "../data/data_raw.h5ad"
ADATA = sc.read(PATH_H5AD)
ADATA.layers["counts"] = ADATA.X.copy()
ADATA_BATCH = ADATA.copy()
ADATA_BATCH.obs["batch"] = 1350 * ["a"] + 1350 * ["b"]

gin.parse_config_file("test-config.gin")


@pytest.mark.parametrize("test_data_path", [(PATH_H5AD), (STR_H5AD)])
def test_load_data(test_data_path):
    assert load_data(test_data_path)


def test_quality_control():
    quality_control(adata=ADATA)


def test_normalization():
    normalization(adata=ADATA)


def test_feature_selection():
    feature_selection(adata=ADATA)


def test_scaling():
    scaling(adata=ADATA)


def test_pca():
    pca(adata=ADATA)


def test_batch_effect_correction():
    pca(ADATA_BATCH, use_highly_variable=False)
    batch_effect_correction(adata=ADATA_BATCH)
