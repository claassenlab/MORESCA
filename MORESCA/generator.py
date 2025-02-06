import argparse
import re
import subprocess
import sys
import warnings
from datetime import datetime
from pathlib import Path

import yaml
from yaml.loader import SafeLoader


def generator(h5ad_path, yaml_path, file_name):
    try:
        with open(yaml_path, "r") as f:
            param_dict = list(yaml.load_all(f, Loader=SafeLoader))[0]
    except FileNotFoundError:
        sys.exit(f"Parameter YAML file {yaml_path} not found.")

    qc_dict = param_dict["QC"]

    doublet_str = ""

    if qc_dict["doublet_removal"]:
        doublet_str = """
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

adata = adata[(~adata.obs.doublet)]"""

    gene_count_filter_str = ""
    match gene_count_flt := qc_dict["n_genes_by_counts"]:
        # Todo: Implement automatic selection of threshold.
        case "auto":
            sys.exit("Auto selection of n_genes_by_counts not implemented.")
        case None:
            pass
        case gene_count_flt if isinstance(gene_count_flt, float):
            gene_count_filter_str = f"""adata[adata.obs.n_genes_by_counts < {gene_count_flt}, :]"
            )"""
        case _:
            sys.exit("Invalid value for n_genes_by_counts.")

    mt_threshold_str = ""

    match qc_dict["mt_threshold"]:
        case mt_value if isinstance(mt_value, (int, float)) and not isinstance(
            mt_value, bool
        ):
            mt_threshold_str = (
                f"adata[adata.obs['pct_counts_mt'] < {qc_dict['mt_threshold']}, :]"
            )
        case "auto":
            sys.exit("Auto selection for mt_threshold not implemented.")
        case None | False:
            pass

    rb_threshold_str = ""

    match qc_dict["rb_threshold"]:
        case rb_value if isinstance(rb_value, (int, float)) and not isinstance(
            rb_value, bool
        ):
            rb_threshold_str = (
                f"adata[adata.obs['pct_counts_ribo'] > {qc_dict['rb_threshold']}, :]"
            )
        case "auto":
            sys.exit("Auto selection for rb_threshold not implemented.")
        case None | False:
            pass

    hb_threshold_str = ""

    match qc_dict["hb_threshold"]:
        case hb_value if isinstance(hb_value, (int, float)) and not isinstance(
            hb_value, bool
        ):
            hb_threshold_str = (
                f"adata[adata.obs['pct_counts_hb'] < {qc_dict['hb_threshold']}, :]"
            )
        case "auto":
            sys.exit("Auto selection for hb_threshold not implemented.")
        case None | False:
            pass

    # sc.pl.highest_expr_genes(adata, n_top=20, show=False, save=True)

    filter_cells_str = f"sc.pp.filter_cells(adata, min_genes={qc_dict['min_genes']})"
    filter_genes_str = f"sc.pp.filter_genes(adata, min_cells={qc_dict['min_cells']})"

    match qc_dict["remove_mt"]:
        case True:
            mito_gene_rm_str = "gene_stack_lst.append(mito_genes)"
        case False | None:
            mito_gene_rm_str = ""
        case _:
            sys.exit("Invalid choice for remove_mt.")

    match qc_dict["remove_rb"]:
        case True:
            ribo_gene_rm_str = "gene_stack_lst.append(ribo_genes)"
        case False | None:
            ribo_gene_rm_str = ""
        case _:
            sys.exit("Invalid choice for remove_rb.")

    match qc_dict["remove_hb"]:
        case True:
            hb_gene_rm_str = "gene_stack_lst.append(hb_genes)"
        case False | None:
            hb_gene_rm_str = ""
        case _:
            sys.exit("Invalid choice for remove_hb.")

    if qc_dict["remove_custom_genes"] is not None:
        warnings.warn(
            "Removing custom genes is not implemented yet. Continue without doing this.",
            category=RuntimeWarning,
        )

    norm_method_str = ""

    match norm_method := qc_dict["normalization"]:
        case "log1pCP10k":
            norm_method_str = """
sc.pp.normalize_total(adata, target_sum=10e4)
sc.pp.log1p(adata)"""
        case "log1PF":
            norm_method_str = """
sc.pp.normalize_total(adata, target_sum=None)
sc.pp.log1p(adata)"""
        case "PFlog1pPF":
            norm_method_str = """
sc.pp.normalize_total(adata, target_sum=None)
sc.pp.log1p(adata)
sc.pp.normalize_total(adata, target_sum=None)"""
        case "pearson_residuals":
            norm_method_str = "sc.experimental.pp.normalize_pearson_residuals(adata)"
        case None | False:
            pass
        case _:
            sys.exit(f"Normalization method {norm_method} not available.")

    feature_number = qc_dict["number_features"]

    feature_method_str = ""

    match feature_method := qc_dict["feature_selection"]:
        case "seurat":
            feature_method_str = (
                f"""sc.pp.highly_variable_genes(adata, flavor="{feature_method}")"""
            )
        case "seurat_v3":
            feature_method_str = f"""sc.pp.highly_variable_genes(adata, flavor="{feature_method}",
            n_top_genes={feature_number}, layer="counts")"""
        case "pearson_residuals":
            feature_method_str = f"""sc.experimental.pp.highly_variable_genes(adata, flavor='{feature_method}',
            n_top_genes={feature_number})"""
        case "anti_correlation":
            feature_method_str = """
warnings.warn(
    'This feature selection is currently only implemented for human data!',
    category=RuntimeWarning,
)
# This is experimental and has to be tested and discussed!
# Todo: Implement mapping for species according to provided YAML.
from anticor_features.anticor_features import get_anti_cor_genes

anti_cor_table = get_anti_cor_genes(
    adata.X.T, adata.var.index.tolist(), species='hsapiens'
)
anti_cor_table.fillna(value=False, axis=None, inplace=True)
adata.var['highly_variable'] = anti_cor_table.selected.copy()"""
        case None | False:
            pass
        case _:
            sys.exit(
                f"Selected feature selection method {feature_method} not available."
            )

    scale_str = "sc.pp.scale(adata)" if qc_dict["scale"] else ""
    batch_dict = param_dict["BatchEffectCorrection"]

    batch_correct_str = ""

    match batch_dict["method"]:
        case "harmony":
            batch_correct_str = """sce.pp.harmony_integrate(
                adata=adata,
                key=batch_dict["batch_key"],
                basis="X_pca",
                adjusted_basis="X_pca",
                max_iter_harmony=50,
            )"""
        case False | None:
            pass
        case _:
            sys.exit("Invalid choice for batch effect correction method.")

    """
    # Preprocessing

    neighbor_dict = param_dict['NeighborhoodGraph']

    # Make this depending on integration choice.
    sc.pp.neighbors(
        adata,
        n_neighbors=neighbor_dict['n_neighbors'],
        n_pcs=neighbor_dict['n_pcs'],
        use_rep='X_pca',
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
        case None:
            print("No clustering done. Exiting.")
            sys.exit(0)
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

    adata.write(Path(RESULT_PATH, "adata_processed.h5ad"))
    """
    dt_string = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    generated_code = f"""# Automatically generated code at {dt_string}
import doubletdetection
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scipy.stats as ss

from anndata import AnnData

def is_outlier(adata: AnnData, metric: str, nmads: int) -> pd.Series(dtype=bool):
    M = adata.obs[metric]
    MAD = ss.median_abs_deviation(M)
    outlier = (M < np.median(M) - nmads * MAD) | (np.median(M) + nmads * MAD < M)
    return outlier

adata = sc.read("{h5ad_path}")
{doublet_str}
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
{gene_count_filter_str}
{mt_threshold_str}
{rb_threshold_str}
{hb_threshold_str}
{gene_count_filter_str}
{filter_cells_str}
{filter_genes_str}

gene_stack_lst = []
gene_stack_lst.append(np.zeros_like(a=adata.var_names))
{mito_gene_rm_str}
{ribo_gene_rm_str}
{hb_gene_rm_str}
remove = np.stack(gene_stack_lst).sum(axis=0).astype(bool)
keep = np.invert(remove)
adata = adata[:, keep]
adata.layers["counts"] = adata.X.copy()
{norm_method_str}
adata.raw = adata
{feature_method_str}
{scale_str}
sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
{batch_correct_str}
"""
    if file_name.suffix != ".py":
        file_name = f"{file_name}.py"

    generated_code = re.sub(r"\n\s*\n", "\n\n", generated_code)

    with open(file_name, mode="w") as file:
        file.write(generated_code)

    # Path to the file you want to format
    file_path = "generated_code.py"

    # Command to run Black on the file
    black_cmd = f"conda run -n scTemplate black {file_path}"

    # Run the command using subprocess
    subprocess.run(black_cmd, shell=True)


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
        "-n",
        "--name",
        type=Path,
        default=Path("generated_code.py"),
        help="Name of the generated Python file.",
    )
    args = parser.parse_args()
    generator(h5ad_path=args.data, yaml_path=args.parameters, file_name=args.name)
