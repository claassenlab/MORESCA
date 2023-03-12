import argparse
import doubletdetection
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scipy.stats as ss
import sys
import warnings
import yaml

from anndata import AnnData
from pathlib import Path
from yaml.loader import SafeLoader


def is_outlier(adata: AnnData, metric: str, nmads: int) -> pd.Series(dtype=bool):
    M = adata.obs[metric]
    MAD = ss.median_abs_deviation(M)
    outlier = (M < np.median(M) - nmads * MAD) | (np.median(M) + nmads * MAD < M)
    return outlier


def run_analysis(h5adPath: Path, yamlPath: Path, figures: bool, verbose: bool) -> None:
    FIGURE_PATH = Path("figures")
    FIGURE_PATH_PRE = Path(FIGURE_PATH, "preQC")
    FIGURE_PATH_POST = Path(FIGURE_PATH, "postQC")
    RESULT_PATH = Path("results")

    FIGURE_PATH.mkdir(exist_ok=True)
    RESULT_PATH.mkdir(exist_ok=True)
    FIGURE_PATH_PRE.mkdir(exist_ok=True)
    FIGURE_PATH_POST.mkdir(exist_ok=True)

    sc.settings.figdir = Path(FIGURE_PATH_PRE)

    try:
        with open(yamlPath, "r") as f:
            param_dict = list(yaml.load_all(f, Loader=SafeLoader))[0]
    except FileNotFoundError:
        sys.exit(f"Parameter YAML file {yamlPath} not found.")

    # Load data.
    adata = sc.read(h5adPath)

    # Quality control
    qc_dict = param_dict["QC"]

    # Todo: This should be aware of the batch key.
    if qc_dict["doublet_removal"]:
        clf = doubletdetection.BoostClassifier(
            n_iters=10,
            clustering_algorithm="phenograph",
            standard_scaling=True,
            pseudocount=0.1,
            n_jobs=-1,
        )

        adata.obs["doublet"] = clf.fit(adata.X).predict(
            p_thresh=1e-16, voter_thresh=0.5
        )
        adata.obs["doublet"] = adata.obs["doublet"].astype(bool)
        adata.obs["doublet_score"] = clf.doublet_score()

        adata = adata[(~adata.obs.doublet)]

    print(f"Cells in raw data: {adata.n_obs}")

    # Quality control - calculate QC covariates
    adata.obs["n_counts"] = adata.X.sum(1)
    adata.obs["log_counts"] = np.log(adata.obs["n_counts"])
    adata.obs["n_genes"] = (adata.X > 0).sum(1)

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.contains(("^RP[SL]"))
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=[20],
        log1p=True,
        inplace=True,
    )

    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )

    # Idea: Make one Figure before QC using matplotlib.pyplot.subplot_mosaic?
    if figures:
        sc.pl.violin(
            adata,
            keys=[
                "log1p_total_counts",
                "log1p_n_genes_by_counts",
                "pct_counts_in_top_20_genes",
            ],
            multi_panel=True,
            show=False,
            save="qc_before_outlier_removal",
        )

    adata = adata[(~adata.obs.outlier)]

    if figures:
        sc.pl.violin(
            adata,
            [
                "n_genes_by_counts",
                "log_counts",
                "pct_counts_mt",
                "pct_counts_ribo",
                "pct_counts_hb",
            ],
            show=False,
            save="preQC",
            multi_panel=True
        )

        sc.pl.scatter(
            adata, x="total_counts", y="n_genes_by_counts", show=False, save="scatter"
        )

    match gene_count_flt := qc_dict["n_genes_by_counts"]:
        # Todo: Implement automatic selection of threshold.
        case "auto" | None:
            pass
        case gene_count_flt if isinstance(gene_count_flt, int):
            adata = adata[adata.obs.n_genes_by_counts < gene_count_flt, :]
        case _:
            sys.exit("Invalid value for n_genes_by_counts.")

    match qc_dict["mt_threshold"]:
        case mt_value if isinstance(mt_value, (int, float)) and not isinstance(
            mt_value, bool
        ):
            adata = adata[adata.obs["pct_counts_mt"] < qc_dict["mt_threshold"], :]
        case "auto" | None:
            # Todo: Implement automatic selection of threshold.
            pass
        case False:
            pass

    match qc_dict["rb_threshold"]:
        case rb_value if isinstance(rb_value, (int, float)) and not isinstance(
            rb_value, bool
        ):
            adata = adata[adata.obs["pct_counts_ribo"] > qc_dict["rb_threshold"], :]
        case "auto" | None:
            # Todo: Implement automatic selection of threshold.
            pass
        case False:
            pass

    match qc_dict["hb_threshold"]:
        case hb_value if isinstance(hb_value, (int, float)) and not isinstance(
            hb_value, bool
        ):
            adata = adata[adata.obs["pct_counts_hb"] < qc_dict["hb_threshold"], :]
        case "auto" | None:
            # Todo: Implement automatic selection of threshold.
            pass
        case False:
            pass

    print(f"Total number of cells: {adata.n_obs}")

    sc.pl.highest_expr_genes(adata, n_top=20, show=False, save=True)

    sc.pp.filter_cells(adata, min_genes=qc_dict["min_genes"])
    print(f"Number of cells after gene filter: {adata.n_obs}")
    sc.pp.filter_genes(adata, min_cells=qc_dict["min_cells"])
    print(f"Number of genes after cell filter: {adata.n_vars}")

    mito_genes = adata.var_names.str.startswith("(?i)MT-")
    ribo_genes = adata.var_names.str.contains(("(?i)^RP[SL]"))
    hb_genes = adata.var_names.str.contains("(?i)^HB[^(P)]")

    gene_stack_lst = []

    if qc_dict["remove_mt"]:
        gene_stack_lst.append(mito_genes)
    if qc_dict["remove_rb"]:
        gene_stack_lst.append(ribo_genes)
    if qc_dict["remove_hb"]:
        gene_stack_lst.append(hb_genes)

    if qc_dict["remove_custom_genes"] is not None:
        warnings.warn(
            "Removing custom genes is not implemented yet. Continue without doing this.",
            category=RuntimeWarning,
        )

    remove = np.stack(gene_stack_lst).sum(axis=0).astype(bool)
    keep = np.invert(remove)
    adata = adata[:, keep]

    adata.layers["counts"] = adata.X.copy()

    sc.settings.figdir = Path(FIGURE_PATH_POST)

    match norm_method := qc_dict["normalization"]:
        case "log1pCP10k":
            sc.pp.normalize_total(adata, target_sum=10e4)
            sc.pp.log1p(adata)
        case "log1PF":
            sc.pp.normalize_total(adata, target_sum=None)
            sc.pp.log1p(adata)
        case "PFlog1pPF":
            sc.pp.normalize_total(adata, target_sum=None)
            sc.pp.log1p(adata)
            sc.pp.normalize_total(adata, target_sum=None)
        case "analytical_pearson":
            sc.experimental.pp.normalize_pearson_residuals(adata)
        case None:
            print("No normalization applied.")
        case _:
            sys.exit(f"Normalization method {norm_method} not available.")

    adata.raw = adata

    feature_number = qc_dict["number_features"]

    match feature_method := qc_dict["feature_selection"]:
        case "seurat":
            sc.pp.highly_variable_genes(adata, flavor=feature_method)
        case "seurat_v3":
            sc.pp.highly_variable_genes(
                adata, flavor=feature_method, n_top_genes=feature_number, layer="counts"
            )
        case "analytical_pearson":
            sc.experimental.pp.highly_variable_genes(
                adata, flavor="pearson_residuals", n_top_genes=feature_number
            )
        case "anti_correlation":
            # This is experimental and has to be tested and discussed!
            # Todo: Implement mapping for species according to provided YAML.
            from anticor_features.anticor_features import get_anti_cor_genes

            anti_cor_table = get_anti_cor_genes(
                adata.X.T, adata.var.index.tolist(), species="hsapiens")
            anti_cor_table.fillna(value=False, axis=None, inplace=True)
            adata.var["highly_variable"] = anti_cor_table.selected.copy()
        case None:
            print("No feature selection applied.")
        case _:
            sys.exit(
                f"Selected feature selection method {feature_method} not available."
            )

    if qc_dict["scale"]:
        sc.pp.scale(adata)

    sc.pp.pca(adata, n_comps=50, use_highly_variable=True)

    batch_dict = param_dict["BatchEffectCorrection"]

    match batch_dict["method"]:
        case "harmony":
            sce.pp.harmony_integrate(
                adata=adata,
                key=batch_dict["batch_key"],
                basis="X_pca",
                adjusted_basis="X_pca",
                max_iter_harmony=50,
            )
        case False | None:
            print("No batch effect correction applied.")
        case _:
            sys.exit("Invalid choice for batch effect correction method.")

    # Preprocessing

    pre_dict = param_dict["Preprocessing"]

    # Make this depending on integration choice.
    sc.pp.neighbors(
        adata,
        n_neighbors=pre_dict["NeighborhoodGraph"]["n_neighbors"],
        n_pcs=pre_dict["NeighborhoodGraph"]["n_pcs"],
        use_rep="X_pca",
        random_state=0,
    )

    # Todo: Should "logreg" be the default?
    if (dge_temp := param_dict["DiffGeneExp"]["method"]) is None:
        dge_method = "wilcoxon"
    else:
        dge_method = dge_temp

    cluster_dict = param_dict["Clustering"]

    match cluster_method := cluster_dict["method"]:
        case "leiden":
            resolution = cluster_dict["resolution"]
            sc.tl.leiden(
                adata=adata,
                resolution=resolution,
                key_added=f"leiden_r{resolution}",
                random_state=0,
            )
        case _:
            sys.exit(f"Clustering method {cluster_method} not available.")

    with warnings.catch_warnings():
        warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
        sc.tl.rank_genes_groups(
            adata,
            f"leiden_r{resolution}",
            method=dge_method,
            use_raw=True,
            key_added=f"leiden_r{resolution}_{dge_method}",
        )

        dedf_leiden = sc.get.rank_genes_groups_df(
            adata, group=None, key=f"leiden_r{resolution}_{dge_method}"
        )
        dedf_leiden.drop("pvals", axis=1, inplace=True)
        # Todo: Should we keep all genes, e.g., for later visualization?
        dedf_leiden = dedf_leiden[dedf_leiden["pvals_adj"] < 0.05]

        with pd.ExcelWriter(
            path=f"results/dge_leiden_r{resolution}_{dge_method}.xlsx"
        ) as writer:
            for cluster_id in dedf_leiden.group.unique():
                df_sub_cl = dedf_leiden[dedf_leiden.group == cluster_id].copy()
                df_sub_cl.to_excel(writer, sheet_name=f"c{cluster_id}")

        del df_sub_cl
        del dedf_leiden

    adata.write(Path(RESULT_PATH, "adata_processed.h5ad"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--data",
        type=Path,
        default=Path("data/data_raw.h5ad"),
        help="Path to the H5AD file.",
    )
    parser.add_argument(
        "-p",
        "--parameters",
        type=Path,
        default=Path("parameters.yml"),
        help="Path to the YAML file.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Set verbosity level.",
    )
    parser.add_argument(
        "-f",
        "--figures",
        action="store_true",
        help="Set whether figures will be generated.",
    )
    args = parser.parse_args()

    run_analysis(
        h5adPath=args.data,
        yamlPath=args.parameters,
        figures=args.figures,
        verbose=args.verbose,
    )
