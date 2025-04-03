import pytest
import scanpy as sc

from MORESCA.utils import remove_cells_by_pct_counts, remove_genes

example_data = sc.datasets.pbmc3k()
example_data.var["mt"] = example_data.var_names.str.contains("(?i)^MT-")
example_data.var["rb"] = example_data.var_names.str.contains(("(?i)^RP[SL]"))
example_data.var["hb"] = example_data.var_names.str.contains(("(?i)^HB[^(P)]"))

mt_genes = example_data.var_names.str.contains("(?i)^MT-")
print(mt_genes.sum())
rb_genes = example_data.var_names.str.contains(("(?i)^RP[SL]"))
print(rb_genes.sum())
hb_genes = example_data.var_names.str.contains("(?i)^HB[^(P)]")
print(hb_genes.sum())

rmv_lst = []

sc.pp.calculate_qc_metrics(
    example_data, qc_vars=["mt", "rb", "hb"], percent_top=[20], log1p=True, inplace=True
)


# TODO: This should include more cases.
@pytest.mark.parametrize(
    "adata, genes, threshold",
    [
        (example_data, "mt", "auto"),
        (example_data, "rb", "auto"),
        (example_data, "hb", "auto"),
        (example_data, "mt", 1),
        (example_data, "rb", 1),
        (example_data, "hb", 1),
        (example_data, "mt", 10),
        (example_data, "rb", 10),
        (example_data, "hb", 10),
        (example_data, "mt", 50),
        (example_data, "rb", 50),
        (example_data, "hb", 50),
    ],
)
def test_remove_cells_by_pct_counts(adata, genes, threshold):
    adata = adata.copy()
    adata.uns["MORESCA"] = dict()
    adata.uns["MORESCA"]["quality_control"] = dict()
    remove_cells_by_pct_counts(adata=adata, genes=genes, threshold=threshold)
    if threshold != "auto":
        if genes == "rb":
            assert adata.obs[f"pct_counts_{genes}"].min() > threshold
        else:
            assert adata.obs[f"pct_counts_{genes}"].max() < threshold
    else:
        # TODO: How should we design this test?
        pass


def test_remove_genes_simple():
    lst0 = [1, 0, 1]
    lst1 = [1, 0, 0]
    lst2 = [0, 1, 1]

    rmv_lst = []

    remove_genes(lst0, rmv_lst, True)
    remove_genes(lst1, rmv_lst, True)
    remove_genes(lst2, rmv_lst, True)

    assert rmv_lst == [[1, 0, 1], [1, 0, 0], [0, 1, 1]]


@pytest.mark.parametrize(
    "gene_lst, rmv_lst, gene_key",
    [(mt_genes, rmv_lst, True), (rb_genes, rmv_lst, False)],
)
def test_remove_genes_real(gene_lst, rmv_lst, gene_key):
    remove_genes(gene_lst=gene_lst, rmv_lst=rmv_lst, gene_key=gene_key)
