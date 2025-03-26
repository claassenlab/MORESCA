import argparse
from pathlib import Path
from typing import List, Union

import gin

import MORESCA
from MORESCA.analysis_steps import (
    batch_effect_correction,
    clustering,
    diff_gene_exp,
    feature_selection,
    load_data,
    neighborhood_graph,
    normalization,
    pca,
    plotting,
    quality_control,
    scaling,
    umap,
)


def run_analysis(
    data_path: Union[Path, List[Path]],
    config_path: Path,
    verbose: bool = False,
    result_path: Union[Path, List[Path]] = Path("results"),
) -> None:
    # Check if data_path and result_path are lists or single paths
    if isinstance(data_path, Path):
        data_path = [data_path]
    if isinstance(result_path, Path):
        result_path = [result_path]
    if len(data_path) != len(result_path) and len(result_path) != 1:
        raise ValueError(
            "Incompatible values for `result_path` and `data_path`: "
            "Either specify a single result path or one per input."
        )

    gin.parse_config_file(config_path)
    sample_id = None

    # Run analysis for each data set
    for i, d_path in enumerate(data_path):
        if len(result_path) == 1:
            res_path = result_path[0]
        else:
            res_path = result_path[i]
        res_path.mkdir(parents=True, exist_ok=True)

        if len(data_path) > 1:
            sample_id = f"s{i + 1:02d}"

        # Run analysis steps
        adata = load_data(d_path)
        adata.raw = adata.copy()
        adata.layers["counts"] = adata.X.copy()
        quality_control(adata=adata, sample_id=sample_id)
        normalization(adata=adata)
        feature_selection(adata=adata)
        adata.layers["unscaled"] = adata.X.copy()
        scaling(adata=adata)
        pca(adata)
        batch_effect_correction(adata=adata)
        neighborhood_graph(adata=adata)
        clustering(adata=adata)
        diff_gene_exp(adata=adata, sample_id=sample_id)
        umap(adata=adata)
        plotting(adata=adata, sample_id=sample_id)

        # Add sample ID if several data sets are processed and output is
        #  saved to the same folder
        sample_str = ""
        if len(data_path) > 1 and len(result_path) == 1:
            sample_str = f"_{sample_id}"
        adata.write(Path(res_path, f"data_processed{sample_str}.h5ad"))


def main():
    parser = argparse.ArgumentParser(
        prog="moresca",
        description="Program: MORESCA (MOdular and REproducible Single-"
        "Cell Analysis)\n"
        f"Version: {MORESCA.__version__}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
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

    args = parser.parse_args()

    run_analysis(
        data_path=args.data,
        result_path=args.output,
        config_path=args.parameters,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
