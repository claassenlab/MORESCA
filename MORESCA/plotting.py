from pathlib import Path
import scanpy as sc
from anndata import AnnData
from matplotlib import pyplot as plt
import seaborn as sns


def plot_qc_vars(adata: AnnData, pre_qc: bool, out_dir: Path) -> None:
    # Copy the AnnData object to ensure that it is not modified
    adata = adata.copy()

    # Calculate QC metrics for raw AnnData if plots pre QC are generated
    if pre_qc:
        adata = adata.raw.to_adata()
        adata.var["mt"] = adata.var_names.str.contains("(?i)^MT-")
        adata.var["rb"] = adata.var_names.str.contains("(?i)^RP[SL]")
        adata.var["hb"] = adata.var_names.str.contains("(?i)^HB(?!EGF|S1L|P1).+")
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "rb", "hb"], log1p=True, inplace=True
        )

    # Plot cell level QC metrics
    qc_vars_cells = [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "pct_counts_rb",
        "pct_counts_hb",
    ]
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(9, 6))
    sc.pl.scatter(
        adata, x="total_counts", y="n_genes_by_counts", ax=axs.flat[0], show=False
    )
    for qc_var, ax in zip(qc_vars_cells, axs.flat[1:]):
        sns.violinplot(adata.obs[qc_var], ax=ax)
        sns.stripplot(adata.obs[qc_var], jitter=0.4, s=1, color="black", ax=ax)
    fig.tight_layout()
    fig.savefig(Path(out_dir, f"qc_vars_cells_{'pre' if pre_qc else 'post'}_qc.png"))
    plt.close()

    # Plot gene level QC metrics
    qc_vars_genes = ["n_cells_by_counts", "pct_dropout_by_counts"]
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(6, 3))
    for qc_var, ax in zip(qc_vars_genes, axs.flat):
        sns.violinplot(adata.var[qc_var], ax=ax)
        sns.stripplot(adata.var[qc_var], jitter=0.4, s=1, color="black", ax=ax)
    fig.tight_layout()
    fig.savefig(Path(out_dir, f"qc_vars_genes_{'pre' if pre_qc else 'post'}_qc.png"))
    plt.close()
