import gin
import scanpy as sc

from MORESCA.analysis_steps import (
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

# PATH_H5AD = Path("../data/data_raw.h5ad")
# STR_H5AD = "../data/data_raw.h5ad"
ADATA = sc.datasets.pbmc3k()
# sc.read(PATH_H5AD)
ADATA.layers["counts"] = ADATA.X.copy()
ADATA_BATCH = ADATA.copy()
ADATA_BATCH.obs["batch"] = 1350 * ["a"] + 1350 * ["b"]

gin.parse_config_file("test-config.gin")


# Todo: Adapt the CI tests to automatically download the data.
# @pytest.mark.parametrize("test_data_path", [(PATH_H5AD), (STR_H5AD)])
# def test_load_data(test_data_path):
#    assert load_data(test_data_path)


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


def test_batch_effect_correction():
    pca(ADATA_BATCH, use_highly_variable=False)
    batch_effect_correction(adata=ADATA_BATCH)


def test_neighborhood_graph():
    feature_selection(adata=ADATA)
    pca(adata=ADATA)
    neighborhood_graph(adata=ADATA)


def test_clustering():
    feature_selection(adata=ADATA)
    pca(adata=ADATA)
    neighborhood_graph(adata=ADATA)
    clustering(adata=ADATA)


def test_diff_gene_exp():
    pca(adata=ADATA, use_highly_variable=False)
    neighborhood_graph(adata=ADATA)
    clustering(adata=ADATA)
    diff_gene_exp(adata=ADATA)
