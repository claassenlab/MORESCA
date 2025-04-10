import gin
import scanpy as sc

from MORESCA.pipeline import (
    batch_effect_correction,
    clustering,
    diff_gene_exp,
    feature_selection,
    neighborhood_graph,
    normalization,
    pca,
    quality_control,
    scaling,
)

ADATA = sc.datasets.pbmc3k()
ADATA.layers["counts"] = ADATA.X.copy()
ADATA_BATCH = ADATA.copy()
ADATA_BATCH.obs["batch"] = 1350 * ["a"] + 1350 * ["b"]

gin.parse_config_file("test-config.gin")


def test_quality_control():
    adata = ADATA.copy()
    quality_control(adata=adata)


def test_normalization():
    adata = ADATA.copy()
    normalization(adata=adata)


def test_feature_selection():
    adata = ADATA.copy()
    feature_selection(adata=adata)


def test_scaling():
    adata = ADATA.copy()
    scaling(adata=adata)


def test_pca():
    adata = ADATA.copy()
    feature_selection(adata=adata)
    pca(adata=adata)


def test_batch_effect_correction():
    adata_batch = ADATA_BATCH.copy()
    pca(adata_batch, use_highly_variable=False)
    batch_effect_correction(adata=adata_batch)


def test_neighborhood_graph():
    adata = ADATA.copy()
    feature_selection(adata=adata)
    pca(adata=adata)
    neighborhood_graph(adata=adata)


def test_clustering():
    adata = ADATA.copy()
    feature_selection(adata=adata)
    pca(adata=adata)
    neighborhood_graph(adata=adata)
    clustering(adata=adata)


def test_diff_gene_exp():
    adata = ADATA.copy()
    pca(adata=adata, use_highly_variable=False)
    neighborhood_graph(adata=adata)
    clustering(adata=adata)
    diff_gene_exp(adata=adata)
