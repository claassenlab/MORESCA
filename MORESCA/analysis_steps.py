import inspect
import warnings
from pathlib import Path
from typing import List, Literal, Optional, Tuple, Union

import doubletdetection
import gin
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scipy.stats as ss
from anndata import AnnData
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

from MORESCA.plotting import plot_qc_vars
from MORESCA.utils import (
    choose_representation,
    remove_cells_by_pct_counts,
    remove_genes,
    store_config_params,
)

try:
    from anticor_features.anticor_features import get_anti_cor_genes

    anti_cor_import_error = False
except ImportError:
    anti_cor_import_error = True
    warnings.warn(
        "Could not import anticor_features,\
        install it using 'pip install anticor-features'"
    )


def is_outlier(adata: AnnData, metric: str, nmads: int) -> pd.Series:
    """
    Check if each value in a given metric column of an AnnData object is an outlier.

    Args:
        adata: An AnnData object containing the data.
        metric: The name of the metric column to check for outliers.
        nmads: The number of median absolute deviations (MADs) away from the median to consider a value as an outlier.

    Returns:
        A pandas Series of boolean values indicating whether each value is an outlier or not.
    """
    data = adata.obs[metric]
    med_abs_dev = ss.median_abs_deviation(data)
    return (data < np.median(data) - nmads * med_abs_dev) | (
        np.median(data) + nmads * med_abs_dev < data
    )


# TODO: Should this happen inplace and just mutate a None-adata?
def load_data(data_path) -> AnnData:
    """
    Load data from a specified file path and return an AnnData object.

    Args:
        data_path: The path to the data file.

    Returns:
        An AnnData object containing the loaded data.

    Raises:
        ValueError: If the file format is unknown.

    Note:
        Currently supports loading of '.h5ad', '.loom', and '.h5' file formats.
    """

    if isinstance(data_path, str):
        data_path = Path(data_path)
    if data_path.is_dir():
        # Todo: Implement this for paths.
        pass
    file_extension = data_path.suffix
    match file_extension:
        case ".h5ad":
            adata = sc.read_h5ad(data_path)
        case ".loom":
            adata = sc.read_loom(data_path)
        case ".h5":
            adata = sc.read_10x_h5(data_path)
        case _:
            try:
                adata = sc.read(data_path)
            except ValueError:
                raise ValueError(f"Unknown file format: {file_extension}")
    adata.var_names_make_unique()
    return adata


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
    n_genes_by_counts: Optional[Union[float, int, str, bool]] = None,
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
        n_genes_by_counts: The threshold for the number of genes detected per cell.
        mt_threshold: The threshold for the percentage of counts in mitochondrial genes.
        rb_threshold: The threshold for the percentage of counts in ribosomal genes.
        hb_threshold: The threshold for the percentage of counts in hemoglobin genes.
        figures: The path to the output directory for the quality control plots.
        pre_qc_plots: Whether to generate plots of QC covariates before quality control or not.
        post_qc_plots: Whether to generate plots of QC covariates after quality control or not.
        doublet_removal: Whether to perform doublet removal or not.
        outlier_removal: Whether to remove outliers or not.
        inplace: Whether to perform the quality control steps in-place or return a modified copy of the AnnData object.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.

    Raises:
        ValueError: If an invalid value is provided for `n_genes_by_counts`.

    Todo:
        - Implement doublet removal for different batches.
        - Implement automatic selection of threshold for `n_genes_by_counts`.
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
        adata.obs["outlier"] = (
            is_outlier(adata, "log1p_total_counts", 5)
            | is_outlier(adata, "log1p_n_genes_by_counts", 5)
            | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
        )

        adata._inplace_subset_obs(~adata.obs.outlier)

    match n_genes_by_counts:
        case n_genes_by_counts if isinstance(n_genes_by_counts, float | int):
            adata._inplace_subset_obs(adata.obs.n_genes_by_counts < n_genes_by_counts)
        case "auto":
            raise NotImplementedError("auto-mode is not implemented.")
        case False | None:
            print("No removal based on n_genes_by_counts.")
        case _:
            raise ValueError("Invalid value for n_genes_by_counts.")

    remove_cells_by_pct_counts(adata=adata, genes="mt", threshold=mt_threshold)
    remove_cells_by_pct_counts(adata=adata, genes="rb", threshold=rb_threshold)
    remove_cells_by_pct_counts(adata=adata, genes="hb", threshold=hb_threshold)

    match min_genes:
        case min_genes if isinstance(min_genes, float | int):
            sc.pp.filter_cells(adata, min_genes=min_genes)
        case False | None:
            print("No removal based on min_genes.")
        case _:
            raise ValueError("Invalid value for min_genes.")

    match min_counts:
        case min_counts if isinstance(min_counts, float | int):
            sc.pp.filter_cells(adata, min_counts=min_counts)
        case False | None:
            print("No removal based on min_counts.")
        case _:
            raise ValueError("Invalid value for min_counts.")

    match max_counts:
        case max_counts if isinstance(max_counts, float | int):
            sc.pp.filter_cells(adata, max_counts=max_counts)
        case False | None:
            print("No removal based on max_counts.")
        case _:
            raise ValueError("Invalid value for max_counts.")

    match min_cells:
        case min_cells if isinstance(min_cells, float | int):
            sc.pp.filter_genes(adata, min_cells=min_cells)
        case False | None:
            print("No removal based on min_cells.")
        case _:
            raise ValueError("Invalid value for min_cells.")

    if post_qc_plots:
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


@gin.configurable
def normalization(
    adata: AnnData,
    apply: bool,
    method: Optional[str] = "log1pPF",
    remove_mt: Optional[bool] = False,
    remove_rb: Optional[bool] = False,
    remove_hb: Optional[bool] = False,
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Normalize gene expression data in an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to apply the normalization steps or not.
        method: The normalization method to use. Available options are:
            - "log1pCP10k": Normalize total counts to 10,000 and apply log1p transformation.
            - "log1pPF": Normalize counts per cell to median of total counts and apply log1p transformation.
            - "PFlog1pPF": Normalize counts per cell to median of total counts, apply log1p transformation, and normalize again using the median of total counts.
            - "analytical_pearson": Normalize using analytical Pearson residuals.
        remove_mt: Whether to remove mitochondrial genes or not.
        remove_rb: Whether to remove ribosomal genes or not.
        remove_hb: Whether to remove hemoglobin genes or not.
        inplace: Whether to perform the normalization steps in-place or return a modified copy of the AnnData object.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.

    Raises:
        ValueError: If an invalid normalization method is provided.
    """

    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=normalization.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    match method:
        case "log1pCP10k":
            sc.pp.normalize_total(adata, target_sum=10e4)
            sc.pp.log1p(adata)
        case "log1pPF":
            sc.pp.normalize_total(adata, target_sum=None)
            sc.pp.log1p(adata)
        case "PFlog1pPF":
            sc.pp.normalize_total(adata, target_sum=None)
            sc.pp.log1p(adata)
            sc.pp.normalize_total(adata, target_sum=None)
        case "analytical_pearson":
            sc.experimental.pp.normalize_pearson_residuals(adata)
        case None | False:
            print("No normalization applied.")
            return None
        case _:
            raise ValueError(f"Normalization method {method} not available.")

    mt_genes = adata.var_names.str.contains("(?i)^MT-")
    rb_genes = adata.var_names.str.contains("(?i)^RP[SL]")
    hb_genes = adata.var_names.str.contains("(?i)^HB[^(P)]")

    gene_stack_lst = []

    remove_genes(gene_lst=mt_genes, rmv_lst=gene_stack_lst, gene_key=remove_mt)
    remove_genes(gene_lst=rb_genes, rmv_lst=gene_stack_lst, gene_key=remove_rb)
    remove_genes(gene_lst=hb_genes, rmv_lst=gene_stack_lst, gene_key=remove_hb)

    # Add zero array in case all three selection are not selected.
    gene_stack_lst.append(np.zeros_like(a=adata.var_names))
    remove = np.stack(gene_stack_lst).sum(axis=0).astype(bool)
    keep = np.invert(remove)
    adata = adata[:, keep]

    if not inplace:
        return adata


@gin.configurable
def feature_selection(
    adata: AnnData,
    apply: bool,
    method: Optional[str] = "seurat",
    number_features: Optional[int] = None,
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Perform feature selection on an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to apply the feature selection steps or not.
        method: The feature selection method to use. Available options are:
            - "seurat": Use Seurat's highly variable genes method.
            - "seurat_v3": Use Seurat v3's highly variable genes method.
            - "analytical_pearson": Use analytical Pearson residuals for feature selection.
            - "anti_correlation": Use anti-correlation method for feature selection (currently only implemented for human data).
        number_features: The number of top features to select (only applicable for certain methods).
        inplace: Whether to perform the feature selection steps in-place or return a modified copy of the AnnData object.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.

    Raises:
        ValueError: If an invalid feature selection method is provided.

    Warnings:
        - The "anti_correlation" method is currently only implemented for human data.
        - If the "anti_correlation" method is selected and the `anticor-features` package is not installed, a warning will be raised.

    Todo:
        - Implement mapping for species according to provided YAML for the "anti_correlation" method.
    """

    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=feature_selection.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    match method:
        case "seurat":
            sc.pp.highly_variable_genes(adata, flavor=method)
        case "seurat_v3":
            sc.pp.highly_variable_genes(
                adata, flavor=method, n_top_genes=number_features, layer="counts"
            )
        case "analytical_pearson":
            sc.experimental.pp.highly_variable_genes(
                adata, flavor="pearson_residuals", n_top_genes=number_features
            )
        case "anti_correlation":
            warnings.warn(
                "This feature selection is currently only implemented for human data!",
                category=RuntimeWarning,
            )
            # This is experimental and has to be tested and discussed!
            # Todo: Implement mapping for species according to provided YAML.

            if anti_cor_import_error:
                warnings.warn(
                    "Anti_cor is not available.\
                    Install it using 'pip install anticor-features."
                )

            anti_cor_table = get_anti_cor_genes(
                adata.X.T, adata.var.index.tolist(), species="hsapiens"
            )
            anti_cor_table.fillna(value=False, axis=None, inplace=True)
            adata.var["highly_variable"] = anti_cor_table.selected.copy()
        case False | None:
            # Todo: Should this be a warning?
            print("No feature selection applied.")
            return None
        case _:
            raise ValueError(
                f"Selected feature selection method {method} not available."
            )

    if not inplace:
        return adata


# Todo: This is just a wrapper, to make the usage of config.gin consistent.
# Does this make sense?
@gin.configurable
def scaling(
    adata: AnnData,
    apply: bool,
    max_value: Optional[Union[int, float]] = None,
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Scale the gene expression data in an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to apply the scaling step or not.
        max_value: The maximum value to which the data will be scaled. If None, the data will be scaled to unit variance.
        inplace: Whether to perform the scaling in-place or return a modified copy of the AnnData object.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.
    """

    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=scaling.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    sc.pp.scale(adata, max_value=max_value)

    if not inplace:
        return adata


@gin.configurable
def pca(
    adata: AnnData,
    apply: bool,
    n_comps: int = 50,
    use_highly_variable: bool = True,
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Perform principal component analysis (PCA) on the gene expression data in an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to apply the PCA or not.
        n_comps: The number of principal components to compute.
        use_highly_variable: Whether to use highly variable genes for PCA computation.
        inplace: Whether to perform the PCA in-place or return a modified copy of the AnnData object.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.
    """

    store_config_params(
        adata=adata,
        analysis_step=inspect.currentframe().f_code.co_name,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    if not inplace:
        adata = adata.copy()

    if use_highly_variable:
        X_data = adata[:, adata.var.highly_variable.values].X
    else:
        X_data = adata.X.copy()

    if n_comps == "auto":
        raise NotImplementedError("auto-mode is not implemented.")
    else:
        pca_ = PCA(n_components=n_comps).fit(X_data)
        X_pca = pca_.transform(X_data)

    n_components = pca_.n_components_

    pca_params = {}
    # Todo: This should be dynamic.
    pca_params["params"] = {
        "zero_center": True,
        "use_highly_variable": use_highly_variable,
        # "mask_var": "highly_variable",
    }
    pca_params["variance"] = pca_.explained_variance_
    pca_params["variance_ratio"] = pca_.explained_variance_ratio_

    adata.obsm["X_pca"] = X_pca[..., :n_comps]
    adata.uns["pca"] = pca_params

    # Code taken from
    # https://github.com/scverse/scanpy/blob/79a5a1c323504cf6df1a19f5c6155b2a0628745e/src/scanpy/preprocessing/_pca/__init__.py#L381
    mask_var = None
    if use_highly_variable:
        mask_var = adata.var["highly_variable"].values

    if mask_var is not None:
        adata.varm["PCs"] = np.zeros(shape=(adata.n_vars, n_comps))
        adata.varm["PCs"][mask_var] = pca_.components_.T
    else:
        adata.varm["PCs"] = pca_.components_.T

    # TODO: Save n_components.
    adata.uns["MORESCA"]["pca"]["n_components"] = n_components

    if not inplace:
        return adata


@gin.configurable
def batch_effect_correction(
    adata: AnnData,
    apply: bool,
    method: Optional[str] = "harmony",
    batch_key: str = "batch",
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Perform batch effect correction on an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to apply the batch effect correction or not.
        method: The batch effect correction method to use. Available options are:
            - "harmony": Use the Harmony algorithm for batch effect correction.
        batch_key: The key in `adata.obs` that identifies the batches.
        inplace: Whether to perform the batch effect correction in-place or return a modified copy of the AnnData object.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.

    Raises:
        ValueError: If an invalid batch effect correction method is provided.

    Note:
        - If `batch_key` is None, no batch effect correction will be performed.
    """

    if batch_key is None:
        return None

    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=batch_effect_correction.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        if not inplace:
            return adata
        return None

    match method:
        case "harmony":
            if "X_pca" not in adata.obsm_keys():
                raise KeyError("X_pca not in adata.obsm. Run PCA first.")
            sce.pp.harmony_integrate(
                adata=adata,
                key=batch_key,
                basis="X_pca",
                adjusted_basis="X_pca_corrected",
                max_iter_harmony=50,
            )
        case False | None:
            print("No batch effect correction applied.")
            return None
        case _:
            raise ValueError("Invalid choice for batch effect correction method.")

    if not inplace:
        return adata


@gin.configurable
def neighborhood_graph(
    adata: AnnData,
    apply: bool,
    n_neighbors: int = 15,
    n_pcs: Optional[int] = None,
    metric: str = "cosine",
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Compute the neighborhood graph for an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to compute the neighborhood graph or not.
        n_neighbors: The number of neighbors to consider for each cell.
        n_pcs: The number of principal components to use for the computation.
        metric: The distance metric to use for computing the neighborhood graph.
        inplace: Whether to perform the computation in-place or return a modified copy of the AnnData object.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.
    """

    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=neighborhood_graph.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    # Compute neighbors graph based on corrected PCA if batch integration was performed, otherwise use PCA
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep="X_pca_corrected"
        if "X_pca_corrected" in adata.obsm_keys()
        else "X_pca"
        if "X_pca" in adata.obsm_keys()
        else None,
        metric=metric,
        random_state=0,
    )

    if not inplace:
        return adata


@gin.configurable
def clustering(
    adata: AnnData,
    apply: bool,
    method: str = "leiden",
    resolution: Union[
        float,
        int,
        List[Union[float, int]],
        Tuple[Union[float, int]],
        Literal["auto"],
    ] = 1.0,
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Perform clustering on an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to perform clustering or not.
        method: The clustering method to use. Available options are:
            - "leiden": Use the Leiden algorithm for clustering.
        resolution: The resolution parameter for the clustering method. Can be a single value or a list of values.
        inplace: Whether to perform the clustering in-place or return a modified copy of the AnnData object.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.

    Raises:
        ValueError: If an invalid clustering method is provided or if the resolution parameter has an invalid type.

    Note:
        - The resolution parameter determines the granularity of the clustering. Higher values result in more fine-grained clusters.
    """

    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=clustering.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    match method:
        case "leiden":
            if (
                not isinstance(resolution, (float, int, list, tuple))
                and resolution != "auto"
            ):
                raise ValueError(
                    f"Invalid type for resolution: {type(resolution)}."
                )

            if isinstance(resolution, (float, int)):
                resolutions = [resolution]
            elif resolution == "auto":
                resolutions = [0.25] + list(np.linspace(0.5, 1.5, 11)) + [2.0]
            else:
                resolutions = resolution

            for res in resolutions:
                sc.tl.leiden(
                    adata=adata,
                    resolution=res,
                    key_added=f"leiden_r{res}",
                    random_state=0,
                )
        case False | None:
            return None
        case _:
            raise ValueError(f"Clustering method {method} not available.")

    # Choose best resolution according to silhouette score
    if len(resolutions) > 1:
        neighbors_params = adata.uns["neighbors"]["params"]
        metric = neighbors_params["metric"]
        use_rep = (
            None
            if "use_rep" not in neighbors_params
            else neighbors_params["use_rep"]
        )
        n_pcs = (
            None
            if "n_pcs" not in neighbors_params
            else neighbors_params["n_pcs"]
        )

        # Use the representation used for neighborhood graph computation
        X = choose_representation(adata, use_rep=use_rep, n_pcs=n_pcs)

        scores = np.zeros(len(resolutions))

        for i, res in enumerate(resolutions):
            scores[i] = silhouette_score(
                X, labels=adata.obs[f"leiden_r{res}"], metric=metric
            )

        best_res = resolutions[np.argmax(scores)]
        adata.obs["leiden"] = adata.obs[f"leiden_r{best_res}"]

        adata.uns["MORESCA"]["clustering"]["best_resolution"] = best_res
        adata.uns["MORESCA"]["clustering"]["resolutions"] = resolutions
        adata.uns["MORESCA"]["clustering"]["silhouette_scores"] = scores

    if not inplace:
        return adata


@gin.configurable(denylist=["sample_id"])
def diff_gene_exp(
    adata: AnnData,
    apply: bool,
    method: str = "wilcoxon",
    groupby: str = "leiden_r1.0",
    use_raw: Optional[bool] = False,
    layer: Optional[str] = "counts",
    corr_method: Literal[
        "benjamini-hochberg", "bonferroni"
    ] = "benjamini-hochberg",
    tables: Optional[Union[Path, str]] = Path("results/"),
    inplace: bool = True,
    sample_id: Optional[str] = None,
) -> Optional[AnnData]:
    """
    Perform differential gene expression analysis on an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to perform differential gene expression analysis or not.
        method: The differential gene expression analysis method to use. Available options are:
            - "wilcoxon": Use the Wilcoxon rank-sum test.
            - "t-test": Use the t-test.
            - "logreg": Use logistic regression.
            - "t-test_overestim_var": Use the t-test with overestimated variance.
        groupby: The key in `adata.obs` that identifies the groups for comparison.
        use_raw: Whether to use the raw gene expression data or not.
        layer: The layer in `adata.layers` to use for the differential gene expression analysis.
        corr_method: The method to use for multiple testing correction.
        tables: The path to the output directory for the differential expression tables.
        inplace: Whether to perform the differential gene expression analysis in-place or return a modified copy of the AnnData object.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.

    Note:
        - The result tables are saved as Excel files if `tables` is True.
        - Only genes with adjusted p-values less than 0.05 are included in the result tables.
    """

    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=diff_gene_exp.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    with warnings.catch_warnings():
        warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

        # Todo: Should "logreg" be the default?
        match method:
            case method if method in {
                "wilcoxon",
                "t-test",
                "logreg",
                "t-test_overestim_var",
            }:
                key_added = f"{groupby}_{method}"
                sc.tl.rank_genes_groups(
                    adata=adata,
                    groupby=groupby,
                    method=method,
                    corr_method=corr_method,
                    use_raw=use_raw,
                    key_added=key_added,
                    layer=layer,
                )

                dedf_leiden = sc.get.rank_genes_groups_df(
                    adata=adata, group=None, key=key_added
                )

                dedf_leiden.drop("pvals", axis=1, inplace=True)
                # Todo: Should we keep all genes, e.g., for later visualization?
                dedf_leiden = dedf_leiden[dedf_leiden["pvals_adj"] < 0.05]

                if tables:
                    if isinstance(tables, str):
                        tables = Path(tables)
                    if sample_id:
                        tables = Path(tables) / f"{sample_id}/"
                    tables.mkdir(parents=True, exist_ok=True)
                    with pd.ExcelWriter(
                        path=f"{tables}/dge_{key_added}.xlsx"
                    ) as writer:
                        for cluster_id in dedf_leiden.group.unique():
                            df_sub_cl = dedf_leiden[
                                dedf_leiden.group == cluster_id
                            ].copy()
                            df_sub_cl.to_excel(writer, sheet_name=f"c{cluster_id}")

            case False | None:
                print("No DGE performed.")
                return None

    if not inplace:
        return adata


@gin.configurable
def umap(
    adata: AnnData, apply: bool, inplace: bool = True
) -> Optional[AnnData]:
    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=plotting.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    sc.tl.umap(adata=adata)

    if not inplace:
        return adata


@gin.configurable(denylist=["sample_id"])
def plotting(
    adata: AnnData,
    apply: bool,
    umap: bool = True,
    path: Path = Path("figures"),
    inplace: bool = True,
    sample_id: Optional[str] = None,
) -> Optional[AnnData]:
    # TODO: Check before merging if we changed adata
    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=plotting.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    path = Path(path)

    # Make subfolder if a sample ID is passed (analysis of multiple samples)
    if sample_id:
        path = path / f"{sample_id}/"
    path.mkdir(parents=True, exist_ok=True)

    if umap:
        sc.pl.umap(adata=adata, show=False)
        plt.savefig(Path(path, "umap.png"))
        plt.close()

    if not inplace:
        return adata
