import warnings
from pathlib import Path
from typing import Optional, Union

import doubletdetection
import gin
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scipy.stats as ss
from anndata import AnnData

from MORESCA.utils import remove_cells_by_pct_counts, remove_genes, store_config_params
from MORESCA.plotting import plot_qc_vars

try:
    from anticor_features.anticor_features import get_anti_cor_genes

    anti_cor_import_error = False
except ImportError:
    anti_cor_import_error = True
    warnings.warn(
        "Could not import anticor_features,\
        install it using 'pip install anticor-features'"
    )


def is_outlier(adata: AnnData, metric: str, nmads: int) -> pd.Series(dtype=bool):
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


# Todo: Should this happen inplace and just mutate a None-adata?
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
        Currently supports loading of '.h5ad', '.loom', and '.hdf5' file formats.
    """

    if isinstance(data_path, str):
        data_path = Path(data_path)
    if data_path.is_dir():
        # Todo: Implement this for paths.
        pass
    file_extension = data_path.suffix
    match file_extension:
        case ".h5ad":
            return sc.read(data_path)
        case ".loom":
            return sc.read_loom(data_path)
        case "hdf5":
            return sc.read_10x_h5(data_path)
        case _:
            raise ValueError(f"Unknown file format: {file_extension}")


@gin.configurable
def quality_control(
    adata: AnnData,
    apply: bool,
    doublet_removal: bool,
    outlier_removal: bool,
    min_genes: int,
    min_cells: int,
    n_genes_by_counts: Optional[Union[float, str, bool]],
    mt_threshold: Optional[Union[int, float, str, bool]],
    rb_threshold: Optional[Union[int, float, str, bool]],
    hb_threshold: Optional[Union[int, float, str, bool]],
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Perform quality control on an AnnData object.

    Args:
        adata: An AnnData object to perform quality control on.
        apply: Whether to apply the quality control steps or not.
        doublet_removal: Whether to perform doublet removal or not.
        outlier_removal: Whether to remove outliers or not.
        min_genes: The minimum number of genes required for a cell to pass quality control.
        min_cells: The minimum number of cells required for a gene to pass quality control.
        n_genes_by_counts: The threshold for the number of genes detected per cell.
        mt_threshold: The threshold for the percentage of counts in mitochondrial genes.
        rb_threshold: The threshold for the percentage of counts in ribosomal genes.
        hb_threshold: The threshold for the percentage of counts in hemoglobin genes.
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
            pass
        case False | None:
            print("No removal based on n_genes_by_counts.")
        case _:
            raise ValueError("Invalid value for n_genes_by_counts.")

    remove_cells_by_pct_counts(adata=adata, genes="mt", threshold=mt_threshold)
    remove_cells_by_pct_counts(adata=adata, genes="rb", threshold=rb_threshold)
    remove_cells_by_pct_counts(adata=adata, genes="hb", threshold=hb_threshold)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    if not inplace:
        return adata


@gin.configurable
def normalization(
    adata: AnnData,
    apply: bool,
    method: str,
    remove_mt: Optional[bool],
    remove_rb: Optional[bool],
    remove_hb: Optional[bool],
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Normalize gene expression data in an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to apply the normalization steps or not.
        method: The normalization method to use. Available options are:
            - "log1pCP10k": Normalize total counts to 10,000 and apply log1p transformation.
            - "log1PF": Normalize counts per cell to median of total counts and apply log1p transformation.
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
        case "log1PF":
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
    method: str,
    number_features=None,
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
    max_value: Optional[Union[int, float]],
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
    use_highly_variable: int = True,
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

    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=pca.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None
    sc.pp.pca(adata, n_comps=n_comps, use_highly_variable=use_highly_variable)

    if not inplace:
        return adata


@gin.configurable
def batch_effect_correction(
    adata: AnnData,
    apply: bool,
    method: str,
    batch_key: str,
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
        return None

    match method:
        case "harmony":
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
    n_neighbors: int,
    n_pcs: int,
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
        else "X_pca",
        metric=metric,
        random_state=0,
    )

    if not inplace:
        return adata


@gin.configurable
def clustering(
    adata: AnnData,
    apply: bool,
    method: str,
    resolution=None,
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
            if not isinstance(resolution, (float, int, list, tuple)):
                raise ValueError(f"Invalid type for resolution: {type(resolution)}.")

            resolutions = (
                [resolution] if isinstance(resolution, (float, int)) else resolution
            )
            for res in resolutions:
                sc.tl.leiden(
                    adata=adata,
                    resolution=res,
                    key_added=f"leiden_r{res}",
                    random_state=0,
                )
        case False | None:
            print("No clustering done. Exiting.")
            return None
        case _:
            raise ValueError(f"Clustering method {method} not available.")

    if not inplace:
        return adata


@gin.configurable
def diff_gene_exp(
    adata: AnnData,
    apply: bool,
    method: str,
    groupby: str,
    use_raw: bool = True,
    layer: str = "counts",
    corr_method: str = "benjamini-hochberg",
    tables: bool = True,
    inplace: bool = True,
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
        tables: Whether to generate result tables or not.
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
                    use_raw=True,
                    key_added=key_added,
                )

                dedf_leiden = sc.get.rank_genes_groups_df(
                    adata=adata, group=None, key=key_added
                )

                dedf_leiden.drop("pvals", axis=1, inplace=True)
                # Todo: Should we keep all genes, e.g., for later visualization?
                dedf_leiden = dedf_leiden[dedf_leiden["pvals_adj"] < 0.05]

                if tables:
                    with pd.ExcelWriter(
                        path=f"results/dge_leiden_r{key_added}.xlsx"
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
    adata: AnnData,
    apply: bool,
    inplace: bool = True,
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


@gin.configurable
def plotting(
    adata: AnnData,
    apply: bool,
    umap: bool,
    pre_qc: bool,
    post_qc: bool,
    path: Path,
    inplace: bool = True,
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
    path.mkdir(exist_ok=True)

    if umap:
        sc.pl.umap(adata=adata, show=False)
        plt.savefig(Path(path, "umap.png"))
        plt.close()

    if pre_qc:
        plot_qc_vars(adata=adata, pre_qc=True, out_dir=path)

    if post_qc:
        plot_qc_vars(adata=adata, pre_qc=False, out_dir=path)

    if not inplace:
        return adata
