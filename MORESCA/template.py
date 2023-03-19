import argparse
import sys
import warnings
from pathlib import Path

import doubletdetection
import gin
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scipy.stats as ss
from analysis_steps import (batch_effect_correction, clustering, diff_gene_exp,
                            feature_selection, load_data, neighborhood_graph,
                            normalization, pca, quality_control, scaling)
from anndata import AnnData
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

    gin.parse_config_file("config.gin")

    adata = load_data(data_path)
    adata.layers["counts"] = adata.X.copy()
    quality_control(adata=adata)
    normalization(adata=adata)
    feature_selection(adata=adata)
    adata.layers["unscaled"] = adata.X.copy()
    adata.raw = adata.copy()
    scaling(adata=adata)
    pca(adata)
    batch_effect_correction(adata=adata)
    neighborhood_graph(adata=adata)
    clustering(adata=adata)
    diff_gene_exp(adata=adata)

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
