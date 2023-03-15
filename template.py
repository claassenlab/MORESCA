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

from typing import Optional

from analysis_steps import load_data, quality_control, normalization, feature_selection, batch_effect_correction, \
    neighborhood_graph, clustering, diff_gene_exp
from utils import remove_cells_by_pct_counts
from utils import remove_genes

from anndata import AnnData
from pathlib import Path
from yaml.loader import SafeLoader


def run_analysis(
    data_path: Path, yaml_path: Path, figures: bool, verbose: bool
) -> None:
    FIGURE_PATH = Path("figures")
    FIGURE_PATH_PRE = Path(FIGURE_PATH, "preQC")
    FIGURE_PATH_POST = Path(FIGURE_PATH, "postQC")
    RESULT_PATH = Path("results")

    FIGURE_PATH.mkdir(exist_ok=True)
    RESULT_PATH.mkdir(exist_ok=True)
    FIGURE_PATH_PRE.mkdir(exist_ok=True)
    FIGURE_PATH_POST.mkdir(exist_ok=True)

    try:
        with open(yaml_path, "r") as f:
            param_dict = list(yaml.load_all(f, Loader=SafeLoader))[0]
    except FileNotFoundError:
        sys.exit(f"Parameter YAML file {yaml_path} not found.")

    qc_dict = param_dict["QC"]
    norm_dict = param_dict["Normalization"]
    feature_dict = param_dict["FeatureSelection"]
    batch_correct_dict = param_dict["BatchEffectCorrection"]
    neighbor_dict = param_dict["NeighborhoodGraph"]
    cluster_dict = param_dict["Clustering"]
    dge_dict = param_dict["DiffGeneExp"]

    adata = load_data(data_path)
    adata.layers["counts"] = adata.X.copy()
    quality_control(adata=adata, **qc_dict)
    normalization(adata=adata, **norm_dict)
    feature_selection(adata=adata, **feature_dict)

    if param_dict["Scale"]:
        # Todo: Better naming.
        adata.layers["unscaled"] = adata.X.copy()
        sc.pp.scale(adata)
    adata.raw = adata
    if param_dict["PCA"]:
        sc.pp.pca(adata, n_comps=50, use_highly_variable=True)

    batch_effect_correction(adata=adata, **batch_correct_dict)
    neighborhood_graph(adata=adata, **neighbor_dict)
    clustering(adata=adata, **cluster_dict)
    diff_gene_exp(adata=adata, **dge_dict)

    # adata.write(Path(RESULT_PATH, "adata_processed.h5ad"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--data",
        type=Path,
        nargs="+",
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
        data_path=args.data,
        yaml_path=args.parameters,
        figures=args.figures,
        verbose=args.verbose,
    )
