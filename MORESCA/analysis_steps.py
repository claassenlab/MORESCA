import warnings
from pathlib import Path
from typing import Optional, Union

import gin
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scipy.stats as ss
from anndata import AnnData

from MORESCA.utils import remove_cells_by_pct_counts, remove_genes


def is_outlier(adata: AnnData, metric: str, nmads: int) -> pd.Series(dtype=bool):
    M = adata.obs[metric]
    MAD = ss.median_abs_deviation(M)
    return (M < np.median(M) - nmads * MAD) | (np.median(M) + nmads * MAD < M)


def load_data(data_path) -> AnnData:
    if isinstance(data_path, str):
        data_path = Path(data_path)
    if data_path.is_dir():
        # Todo: Implement this for paths.
        pass
    file_extension = data_path.suffix
    if file_extension == ".h5ad":
        return sc.read(data_path)
    elif file_extension == ".loom":
        return sc.read_loom(data_path)
    elif file_extension == ".hdf5":
        return sc.read_10x_h5(data_path)
    else:
        raise ValueError(f"Unknown file format: {file_extension}")


@gin.configurable
def quality_control(
    adata: AnnData,
    doublet_removal: bool,
    outlier_removal: bool,
    min_genes: int,
    min_cells: int,
    n_genes_by_counts: Optional[Union[float, str, bool]],
    mt_threshold: Optional[Union[int, float, str, bool]],
    rb_threshold: Optional[Union[int, float, str, bool]],
    hb_threshold: Optional[Union[int, float, str, bool]],
    remove_mt: Optional[bool],
    remove_rb: Optional[bool],
    remove_hb: Optional[bool],
    # Todo: Only for testing, this should be a list.
    remove_custom_genes: Optional[bool],
    inplace: bool = True,
    save: Union[Path, str, bool] = False,
) -> Optional[AnnData]:
    if not inplace:
        adata = adata.copy()

    # Todo: This should be aware of the batch key.
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

        adata = adata[(~adata.obs.doublet)]

    # Quality control - calculate QC covariates
    adata.obs["n_counts"] = adata.X.sum(1)
    adata.obs["log_counts"] = np.log(adata.obs["n_counts"])
    adata.obs["n_genes"] = (adata.X > 0).sum(1)

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["rb"] = adata.var_names.str.contains(("^RP[SL]"))
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "rb", "hb"],
        percent_top=[20],
        log1p=True,
        inplace=True,
    )

    if outlier_removal:
        adata.obs["outlier"] = (
            is_outlier(adata, "log1p_total_counts", 5)
            | is_outlier(adata, "log1p_n_genes_by_counts", 5)
            | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
        )

        adata = adata[(~adata.obs.outlier)]

    match n_genes_by_counts:
        case n_genes_by_counts if isinstance(n_genes_by_counts, float):
            adata = adata[adata.obs.n_genes_by_counts < n_genes_by_counts, :]
        # Todo: Implement automatic selection of threshold.
        case "auto":
            pass
        case False | None:
            print("No removal based on n_genes_by_counts.")
        case _:
            ValueError("Invalid value for n_genes_by_counts.")

    # Todo: Test this!
    remove_cells_by_pct_counts(adata=adata, genes="mt", threshold=mt_threshold)
    remove_cells_by_pct_counts(adata=adata, genes="rb", threshold=rb_threshold)
    remove_cells_by_pct_counts(adata=adata, genes="hb", threshold=hb_threshold)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    mt_genes = adata.var_names.str.startswith("(?i)MT-")
    rb_genes = adata.var_names.str.contains(("(?i)^RP[SL]"))
    hb_genes = adata.var_names.str.contains("(?i)^HB[^(P)]")

    gene_stack_lst = []

    remove_genes(gene_lst=mt_genes, rmv_lst=gene_stack_lst, gene_key=remove_mt)
    remove_genes(gene_lst=rb_genes, rmv_lst=gene_stack_lst, gene_key=remove_rb)
    remove_genes(gene_lst=hb_genes, rmv_lst=gene_stack_lst, gene_key=remove_hb)

    if (remove_custom_genes is not None) or remove_custom_genes:
        warnings.warn(
            "Removing custom genes is not implemented yet. Continue without doing this.",
            category=RuntimeWarning,
        )

    # Add zero array in case all three selection are not selected.
    gene_stack_lst.append(np.zeros_like(a=adata.var_names))
    remove = np.stack(gene_stack_lst).sum(axis=0).astype(bool)
    keep = np.invert(remove)
    adata = adata[:, keep]

    if save:
        if isinstance(save, Path | str):
            adata.write(save)
        else:
            adata.write("results/post_qc.h5ad")

    if not inplace:
        return adata


@gin.configurable
def normalization(
    adata: AnnData, method: str, inplace: bool = True, save: bool = False
) -> Optional[AnnData]:
    if not inplace:
        adata = adata.copy()

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
            ValueError(f"Normalization method {method} not available.")

    if save:
        if isinstance(save, Path | str):
            adata.write(save)
        else:
            adata.write("results/normalized.h5ad")

    if not inplace:
        return adata


@gin.configurable
def feature_selection(
    adata: AnnData,
    method: str,
    number_features=None,
    inplace: bool = True,
    save: bool = False,
) -> Optional[AnnData]:
    if not inplace:
        adata = adata.copy()

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
            from anticor_features.anticor_features import get_anti_cor_genes

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
            ValueError(f"Selected feature selection method {method} not available.")

    if save:
        if isinstance(save, Path | str):
            adata.write(save)
        else:
            # Todo: Do we have to save here?
            adata.write("results/feature_selection.h5ad")

    if not inplace:
        return adata


# Todo: This is just a wrapper, to make the usage of config.gin consistent Does this make sense?
@gin.configurable
def scaling(
    adata: AnnData,
    apply: bool,
    max_value: Optional[Union[int, float]],
    inplace: bool = True,
    save: bool = False,
)-> Optional[AnnData]:
    if not inplace:
        adata = adata.copy()

    if not apply:
        return None

    sc.pp.scale(adata, max_value=max_value)

    if save:
        if isinstance(save, Path | str):
            adata.write(save)
        else:
            # Todo: Do we have to save here?
            adata.write("results/scaled.h5ad")

    if not inplace:
        return adata


@gin.configurable
def pca(
    adata: AnnData,
    apply: bool,
    n_comps: int = 50,
    use_highly_variable: int = True,
    inplace: bool = True,
    save: bool = False,
)-> Optional[AnnData]:
    if not inplace:
        adata = adata.copy()

    if not apply:
        return None
    sc.pp.pca(adata, n_comps=n_comps, use_highly_variable=use_highly_variable)

    if save:
        if isinstance(save, Path | str):
            adata.write(save)
        else:
            # Todo: Do we have to save here?
            adata.write("results/pca.h5ad")

    if not inplace:
        return adata


@gin.configurable
def batch_effect_correction(
    adata: AnnData,
    method: str,
    batch_key: str,
    inplace: bool = True,
    save: bool = False,
) -> Optional[AnnData]:
    if batch_key is None:
        return None

    if not inplace:
        adata = adata.copy()

    match method:
        case "harmony":
            sce.pp.harmony_integrate(
                adata=adata,
                key=batch_key,
                basis="X_pca",
                # Todo: Should this be a different layer?
                adjusted_basis="X_pca",
                max_iter_harmony=50,
            )
        case False | None:
            print("No batch effect correction applied.")
            return None
        case _:
            ValueError("Invalid choice for batch effect correction method.")

    if save:
        if isinstance(save, Path | str):
            adata.write(save)
        else:
            adata.write("results/batch_corrected.h5ad")

    if not inplace:
        return adata


@gin.configurable
def neighborhood_graph(
    adata: AnnData,
    n_neighbors: int,
    n_pcs: int,
    inplace: bool = True,
    save: bool = False,
) -> Optional[AnnData]:
    if not inplace:
        adata = adata.copy()

    # Make this depending on integration choice.
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep="X_pca",
        random_state=0,
    )

    if save:
        if isinstance(save, Path | str):
            adata.write(save)
        else:
            adata.write("results/neighborhood_graph.h5ad")

    if not inplace:
        return adata


@gin.configurable
def clustering(
    adata: AnnData,
    method: str,
    resolution=None,
    inplace: bool = True,
    save: bool = False,
) -> Optional[AnnData]:
    if not inplace:
        adata = adata.copy()

    match method:
        case "leiden":
            resolution = resolution
            sc.tl.leiden(
                adata=adata,
                resolution=resolution,
                key_added=f"leiden_r{resolution}",
                random_state=0,
            )
        case False | None:
            print("No clustering done. Exiting.")
            return None
        case _:
            ValueError(f"Clustering method {method} not available.")

    if save:
        if isinstance(save, Path | str):
            adata.write(save)
        else:
            adata.write("results/clustered.h5ad")

    if not inplace:
        return adata


@gin.configurable
def diff_gene_exp(
    adata: AnnData,
    method: str,
    groupby: str,
    use_raw: bool = True,
    tables: bool = True,
    inplace: bool = True,
    save: bool = False,
) -> Optional[AnnData]:
    if not inplace:
        adata = adata.copy()

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

    if save:
        if isinstance(save, Path | str):
            adata.write(save)
        else:
            adata.write("results/DGE.h5ad")

    if not inplace:
        return adata
