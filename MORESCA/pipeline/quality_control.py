import logging
from pathlib import Path
from typing import Optional, Union

import doubletdetection
import gin
import numpy as np
import scanpy as sc
from anndata import AnnData

from MORESCA.plotting import plot_qc_vars
from MORESCA.utils import is_outlier, remove_cells_by_pct_counts, store_config_params

log = logging.getLogger(__name__)


@gin.configurable(denylist=["sample_id"])
def quality_control(
    adata: AnnData,
    apply: bool,
    doublet_removal: bool = False,
    outlier_removal: bool = False,
    min_genes: Optional[Union[float, int, bool]] = None,
    min_counts: Optional[Union[float, int, bool]] = None,
    max_counts: Optional[Union[float, int, bool]] = None,
    min_cells: Optional[Union[float, int, bool]] = None,
    max_genes: Optional[Union[float, int, str, bool]] = None,
    mt_threshold: Optional[Union[int, float, str, bool]] = None,
    rb_threshold: Optional[Union[int, float, str, bool]] = None,
    hb_threshold: Optional[Union[int, float, str, bool]] = None,
    figures: Optional[Union[Path, str]] = None,
    pre_qc_plots: Optional[bool] = None,
    post_qc_plots: Optional[bool] = None,
    inplace: bool = True,
    sample_id: Optional[str] = None,
) -> Optional[AnnData]:
    """
    Perform quality control on an AnnData object.

    Args:
        adata: An AnnData object to perform quality control on.
        apply: Whether to apply the quality control steps or not.
        min_genes: The minimum number of genes required for a cell to pass quality control.
        min_counts: The minimum total counts required for a cell to pass quality control.
        max_counts: The maximum total counts allowed for a cell to pass quality control.
        min_cells: The minimum number of cells required for a gene to pass quality control.
        max_genes: The threshold for the number of genes detected per cell.
        mt_threshold: The threshold for the percentage of counts in mitochondrial genes.
        rb_threshold: The threshold for the percentage of counts in ribosomal genes.
        hb_threshold: The threshold for the percentage of counts in hemoglobin genes.
        figures: The path to the output directory for the quality control plots.
        pre_qc_plots: Whether to generate plots of QC covariates before quality control or not.
        post_qc_plots: Whether to generate plots of QC covariates after quality control or not.
        doublet_removal: Whether to perform doublet removal or not.
        outlier_removal: Whether to remove outliers or not.
        inplace: Whether to perform the quality control steps in-place or return a modified copy of the AnnData object.
        sample_id: Sample ID for subfolder creation in the output directory. Not gin-configurable.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.

    Raises:
        ValueError: If an invalid value is provided for `max_genes`, `min_genes`, `min_counts`, `max_counts`, or `min_cells`.

    Todo:
        - Implement doublet removal for different batches.
        - Implement automatic selection of threshold for `max_genes`.
    """

    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=quality_control.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    log.info("Performing quality control.")

    # Quality control - calculate QC covariates
    adata.obs["n_counts"] = adata.X.sum(1)
    adata.obs["log_counts"] = np.log(adata.obs["n_counts"])
    adata.obs["n_genes"] = (adata.X > 0).sum(1)

    adata.var["mt"] = adata.var_names.str.contains("(?i)^MT-")
    adata.var["rb"] = adata.var_names.str.contains("(?i)^RP[SL]")
    adata.var["hb"] = adata.var_names.str.contains("(?i)^HB(?!EGF|S1L|P1).+")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "rb", "hb"], percent_top=[20], log1p=True, inplace=True
    )

    if pre_qc_plots:
        log.debug("Creating pre-QC plots.")
        # Make default directory if figures is None or empty string
        if not figures:
            figures = "figures/"
        if isinstance(figures, str):
            figures = Path(figures)

        # Make subfolder if a sample ID is passed (analysis of multiple samples)
        if sample_id:
            figures = figures / f"{sample_id}/"

        figures.mkdir(parents=True, exist_ok=True)
        plot_qc_vars(adata, pre_qc=True, out_dir=figures)

    if doublet_removal:
        log.debug("Performing doublet removal.")
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

        adata._inplace_subset_obs(~adata.obs.doublet)

    if outlier_removal:
        log.debug("Performing outlier removal.")
        adata.obs["outlier"] = (
            is_outlier(adata, "log1p_total_counts", 5)
            | is_outlier(adata, "log1p_n_genes_by_counts", 5)
            | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
        )

        adata._inplace_subset_obs(~adata.obs.outlier)

    match max_genes:
        case max_genes if isinstance(max_genes, float | int):
            log.debug(f"Removing cells with more than {max_genes} genes.")
            sc.pp.filter_cells(adata, max_genes=max_genes)
        case "auto":
            raise NotImplementedError("auto-mode is not implemented.")
        case False | None:
            log.debug("No removal based on max_genes.")
        case _:
            raise ValueError("Invalid value for max_genes.")

    remove_cells_by_pct_counts(adata=adata, genes="mt", threshold=mt_threshold)
    remove_cells_by_pct_counts(adata=adata, genes="rb", threshold=rb_threshold)
    remove_cells_by_pct_counts(adata=adata, genes="hb", threshold=hb_threshold)

    match min_genes:
        case min_genes if isinstance(min_genes, float | int):
            log.debug(f"Removing cells with less than {min_genes} genes.")
            sc.pp.filter_cells(adata, min_genes=min_genes)
        case False | None:
            log.debug("No removal based on min_genes.")
        case _:
            raise ValueError("Invalid value for min_genes.")

    match min_counts:
        case min_counts if isinstance(min_counts, float | int):
            log.debug(f"Removing cells with less than {min_counts} counts.")
            sc.pp.filter_cells(adata, min_counts=min_counts)
        case False | None:
            log.debug("No removal based on min_counts.")
        case _:
            raise ValueError("Invalid value for min_counts.")

    match max_counts:
        case max_counts if isinstance(max_counts, float | int):
            log.debug(f"Removing cells with more than {max_counts} counts.")
            sc.pp.filter_cells(adata, max_counts=max_counts)
        case False | None:
            log.debug("No removal based on max_counts.")
        case _:
            raise ValueError("Invalid value for max_counts.")

    match min_cells:
        case min_cells if isinstance(min_cells, float | int):
            log.debug(f"Removing genes expressed in less than {min_cells} cells.")
            sc.pp.filter_genes(adata, min_cells=min_cells)
        case False | None:
            log.debug("No removal based on min_cells.")
        case _:
            raise ValueError("Invalid value for min_cells.")

    if post_qc_plots:
        log.debug("Creating post-QC plots.")
        # Make default directory if figures is None or empty string
        if not figures:
            figures = "figures/"
        if isinstance(figures, str):
            figures = Path(figures)

        # Make subfolder if a sample ID is passed (analysis of multiple samples)
        if not pre_qc_plots and sample_id:
            figures = figures / f"{sample_id}/"
        figures.mkdir(parents=True, exist_ok=True)
        plot_qc_vars(adata, pre_qc=False, out_dir=figures)

    if not inplace:
        return adata
