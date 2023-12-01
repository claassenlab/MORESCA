import argparse
from pathlib import Path

import gin
from MORESCA.analysis_steps import (
    batch_effect_correction,
    clustering,
    diff_gene_exp,
    feature_selection,
    load_data,
    neighborhood_graph,
    normalization,
    pca,
    quality_control,
    scaling,
    umap,
    plotting,
)


def run_analysis(
    data_path: Path,
    config_path: Path,
    figures: bool,
    verbose: bool,
    result_path: Path = Path("results"),
) -> None:
    result_path.mkdir(exist_ok=True)

    gin.parse_config_file(config_path)

    adata = load_data(data_path)
    adata.raw = adata.copy()
    adata.layers["counts"] = adata.X.copy()
    quality_control(adata=adata)
    normalization(adata=adata)
    feature_selection(adata=adata)
    adata.layers["unscaled"] = adata.X.copy()
    scaling(adata=adata)
    pca(adata)
    batch_effect_correction(adata=adata)
    neighborhood_graph(adata=adata)
    clustering(adata=adata)
    diff_gene_exp(adata=adata)
    umap(adata=adata)
    plotting(adata=adata)

    adata.write(Path(result_path, "data_processed.h5ad"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--data",
        type=Path,
        nargs="+",
        default=Path("data/data_raw.h5ad"),
        help="Path to the data.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        nargs="+",
        default=Path("results"),
        help="Path to the processed output data.",
    )
    parser.add_argument(
        "-p",
        "--parameters",
        type=Path,
        default=Path("config.gin"),
        help="Path to the config.gin.",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Set verbosity level."
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
        result_path=args.output,
        config_path=args.parameters,
        figures=args.figures,
        verbose=args.verbose,
    )
