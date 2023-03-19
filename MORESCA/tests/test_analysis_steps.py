import pytest
import scanpy as sc

from MORESCA import analysis_steps
from pathlib import Path

PATH_H5AD = Path("../data/data_raw.h5ad")
STR_H5AD = "../data/data_raw.h5ad"
ADATA = sc.read(PATH_H5AD)


# Todo: Add cases for each supported format.
@pytest.mark.parametrize("test_data_path", [(PATH_H5AD), (STR_H5AD)])
def test_load_data(test_data_path):
    assert analysis_steps.load_data(test_data_path)


@pytest.mark.parametrize(
    """adata,
    doublet_removal,
    outlier_removal,
    min_genes,
    min_cells,
    n_genes_by_counts,
    mt_threshold,
    rb_threshold,
    hb_threshold,
    remove_mt,
    remove_rb,
    remove_hb,
    remove_custom_genes
    """,
    [
        (
            ADATA,
            False,
            False,
            200,
            10,
            None,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
        )
    ],
)
def test_quality_control(
    adata,
    doublet_removal,
    outlier_removal,
    min_genes,
    min_cells,
    n_genes_by_counts,
    mt_threshold,
    rb_threshold,
    hb_threshold,
    remove_mt,
    remove_rb,
    remove_hb,
    remove_custom_genes,
):
    assert analysis_steps.quality_control()
