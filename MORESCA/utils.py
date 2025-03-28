from typing import Optional, Union

import numpy as np
import scanpy as sc
import scipy.stats as ss
from anndata import AnnData
from scipy.sparse import csr_matrix


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * ss.median_abs_deviation(M)) | (
        np.median(M) + nmads * ss.median_abs_deviation(M) < M
    )
    return outlier


def is_passing_upper(data, nmads, upper_limit=100):
    """
    Check if each value in the given data array is passing the upper limit.

    Args:
        data: An array of values to check.
        nmads: The number of median absolute deviations (MADs) away from the median to define the upper bound.
        upper_limit: The upper limit value. Defaults to 0.

    Returns:
        A boolean array indicating whether each value is passing the upper limit or not.
    """

    med = np.median(data)
    mad = ss.median_abs_deviation(data)
    upper_bound = min(med + nmads * mad, upper_limit)
    return data <= upper_bound


def is_passing_lower(data, nmads, lower_limit=0):
    """
    Check if each value in the given data array is passing the lower limit.

    Args:
        data: An array of values to check.
        nmads: The number of median absolute deviations (MADs) away from the median to define the lower bound.
        lower_limit: The lower limit value.

    Returns:
        A boolean array indicating whether each value is passing the lower limit or not.
    """

    med = np.median(data)
    mad = ss.median_abs_deviation(data)
    lower_bound = max(med - nmads * mad, lower_limit)
    return data >= lower_bound


def remove_cells_by_pct_counts(
    adata: AnnData,
    genes: str,
    threshold: Optional[Union[int, float, str, bool]],
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Remove cells from an AnnData object based on the percentage of counts in specific gene categories.

    Args:
        adata: An AnnData object containing the gene expression data.
        genes: The gene category to filter cells based on. Accepted values are "mt", "rb", and "hb".
        threshold: The threshold value for the percentage of counts. Cells with a percentage above or below this threshold will be removed.
        inplace: Whether to perform the cell removal in-place or return a modified copy of the AnnData object.
        save: Whether to save the modified AnnData object to a file. If a path or file name is provided, the AnnData object will be saved.

    Returns:
        If `inplace` is True and `save` is False, returns None. Otherwise, returns a modified copy of the AnnData object.

    Raises:
        ValueError: If an invalid gene category is provided or if an error occurs during the filtering process.

    Todo:
        - Implement automatic selection of threshold for "auto" option.
    """

    if not inplace:
        adata = adata.copy()

    if genes not in ["mt", "rb", "hb"]:
        raise ValueError(
            f"{genes} is not selectable. Accepted values are ['mt', 'rb'', 'hb']"
        )

    # Should we calculate the auto-thresholds here so they are not affected by other values?
    # E.g., if we filter MT manually, that resulting distribution will be changed, thus
    # the auto-mode results for RB might be differnt compared to running auto for MT and RB.

    adata_aux = adata.copy()

    adata_aux.var["mt"] = adata_aux.var_names.str.contains("(?i)^MT-")
    adata_aux.var["rb"] = adata_aux.var_names.str.contains("(?i)^RP[SL]")
    adata_aux.var["hb"] = adata_aux.var_names.str.contains(
        "(?i)^HB(?!EGF|S1L|P1).+"
    )

    sc.pp.calculate_qc_metrics(
        adata_aux, qc_vars=["mt", "rb", "hb"], percent_top=None, inplace=True
    )

    # TODO: Should we use nmads=3 or nmads=5 here?
    rb_pass = is_passing_lower(adata_aux.obs["pct_counts_rb"], nmads=3)
    hb_pass = is_passing_upper(adata_aux.obs["pct_counts_hb"], nmads=3)

    rb_auto_threshold = adata_aux[rb_pass].obs["pct_counts_rb"].min()
    hb_auto_threshold = adata_aux[hb_pass].obs["pct_counts_hb"].max()

    match threshold:
        case threshold if isinstance(
            threshold, (int, float)
        ) and not isinstance(threshold, bool):
            if genes == "rb":
                adata._inplace_subset_obs(
                    adata.obs[f"pct_counts_{genes}"] > threshold
                )
            else:
                adata._inplace_subset_obs(
                    adata.obs[f"pct_counts_{genes}"] < threshold
                )
        case "auto":
            if genes == "mt":
                ddqc(adata)
            elif genes == "rb":
                adata._inplace_subset_obs(
                    adata.obs[f"pct_counts_{genes}"] > rb_auto_threshold
                )
            elif genes == "hb":
                adata._inplace_subset_obs(
                    adata.obs[f"pct_counts_{genes}"] < hb_auto_threshold
                )
            else:
                raise ValueError(
                    f"Auto selection for {genes}_threshold not implemented."
                )
        case False | None:
            print(f"No {genes} filter applied.")
        case _:
            raise ValueError("Error.")

    if not inplace:
        return adata


# TODO: Is this the best way to do it? Manipulating the list inplace feels like a gotcha.
def remove_genes(gene_lst: list, rmv_lst: list, gene_key) -> None:
    """
    Remove genes from a list based on a specified condition.

    Args:
        gene_lst: A list of genes to check.
        rmv_lst: A list to store the genes to be removed.
        gene_key: A condition to determine whether to remove the genes or not.

    Returns:
        None.

    Raises:
        ValueError: If an invalid choice is provided for `gene_key`.
    """

    match gene_key:
        case True:
            rmv_lst.append(gene_lst)
        case False | None:
            pass
        case _:
            raise ValueError("Invalid choice for gene_key.")


def ddqc(adata: AnnData, inplace: bool = True) -> Optional[AnnData]:
    """
    Perform Data-Driven Quality Control (DDQC) on an AnnData object. Described in the publication

    Biology-inspired data-driven quality control for scientific discovery in single-cell transcriptomics

    https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02820-w

    Args:
        adata: An AnnData object containing the gene expression data.
        inplace: Whether to perform the DDQC in-place or return a modified copy of the AnnData object.
        save: Whether to save the modified AnnData object to a file. If True, the AnnData object will be saved.

    Returns:
        If `inplace` is True and `save` is False, returns None. Otherwise, returns a modified copy of the AnnData object.

    Note:
        - The resulting AnnData object will only contain cells that pass the quality control checks.
        - We modify the published method by only considering the MT threshold.

    Todo:
        - The clustering here is possibly not equivalent to our later clustering routine. Should this match?
    """

    if not inplace:
        adata = adata.copy()

    adata_copy = adata.copy()

    sc.pp.calculate_qc_metrics(
        adata_copy, qc_vars=["mt"], percent_top=None, inplace=True
    )
    adata_copy._inplace_subset_obs(adata_copy.obs.pct_counts_mt <= 80)

    # TODO: can this be removed?
    adata_copy.layers["counts"] = adata_copy.X.copy()
    sc.pp.normalize_total(adata_copy, target_sum=1e4)
    sc.pp.log1p(adata_copy)
    sc.pp.highly_variable_genes(
        adata_copy, flavor="seurat_v3", n_top_genes=2000, layer="counts"
    )
    sc.pp.scale(adata_copy)
    sc.tl.pca(adata_copy)
    sc.pp.neighbors(adata_copy, n_neighbors=20, n_pcs=50, metric="euclidean")
    sc.tl.leiden(adata_copy, resolution=1.4)

    # Directly apply the quality control checks and create the 'passed' mask
    passed = np.ones(adata_copy.n_obs, dtype=bool)

    cellwise_mt_threshold = np.zeros(adata_copy.n_obs, dtype=float)
    # cellwise_n_genes_counts_threshold = np.zeros(adata_copy.n_obs, dtype=float)
    # cellwise_total_counts_threshold = np.zeros(adata_copy.n_obs, dtype=float)

    for cluster in adata_copy.obs["leiden"].unique():
        indices = adata_copy.obs["leiden"] == cluster
        pct_counts_mt_cluster = adata_copy.obs.loc[
            indices, "pct_counts_mt"
        ].values

        passing_mask_mt = is_passing_upper(pct_counts_mt_cluster, nmads=3)
        """
        total_counts_cluster = adata_copy.obs.loc[
            indices, "total_counts"
        ].values
        n_genes_cluster = adata_copy.obs.loc[
            indices, "n_genes_by_counts"
        ].values
        passing_mask_counts = is_passing_lower(
            total_counts_cluster, nmads=3, lower_limit=0
        )
        passing_mask_genes = is_passing_lower(
            n_genes_cluster, nmads=3, lower_limit=200
        )

        total_thresh_ = adata_copy.obs["total_counts"][indices][
            passing_mask_counts
        ].min()
        genes_thresh_ = adata_copy.obs["n_genes_by_counts"][indices][
            passing_mask_genes
        ].min()

        cellwise_n_genes_counts_threshold[indices][passing_mask_counts] = (
            total_thresh_
        )
        cellwise_total_counts_threshold[indices][passing_mask_genes] = (
            genes_thresh_
        )
        """

        mt_thresh_ = adata_copy.obs["pct_counts_mt"][indices][
            passing_mask_mt
        ].max()

        cellwise_mt_threshold[indices][passing_mask_mt] = mt_thresh_

        # & passing_mask_counts & passing_mask_genes
        passed[indices] = passing_mask_mt

    cellwise_mt_threshold = cellwise_mt_threshold[passed]
    # cellwise_n_genes_counts_threshold = cellwise_n_genes_counts_threshold[
    #    passed
    # ]
    # cellwise_total_counts_threshold = cellwise_total_counts_threshold[passed]

    adata.uns["MORESCA"]["quality_control"]["mt_threshold_per_cell"] = (
        cellwise_mt_threshold
    )
    passed = adata_copy[passed].obs_names
    adata._inplace_subset_obs(passed.values)

    if not inplace:
        return adata


def store_config_params(
    adata: AnnData, analysis_step: str, apply: bool, params: dict
) -> None:
    """
    Store configuration parameters for an analysis step in the AnnData object.

    Args:
        adata: An AnnData object.
        analysis_step: The name of the analysis step.
        apply: Whether the analysis step is applied or not.
        params: A dictionary of configuration parameters for the analysis step.

    Returns:
        None.

    Note:
        - The configuration parameters are stored in the `uns` attribute of the AnnData object under the key "MORESCA".
        - The configuration parameters are stored as a dictionary under the analysis step name.
        - If `apply` is False, the configuration parameters are set to False or None, depending on the key.
    """

    # Create a dictionary for storing config parameters
    uns_key = "MORESCA"
    if uns_key not in adata.uns_keys():
        adata.uns[uns_key] = {}
    adata.uns[uns_key][analysis_step] = {}

    # Store config parameters, depending on whether the step is applied or not
    if not apply:
        params = {
            key: (False if key == "apply" else None)
            for key, val in params.items()
        }
    adata.uns[uns_key][analysis_step] = {
        key: (list(val) if isinstance(val, tuple) else val)
        for key, val in params.items()
    }


def choose_representation(
    adata: AnnData,
    use_rep: Optional[str] = None,
    n_pcs: Optional[int] = None,
    **kwargs,
) -> Union[np.ndarray, csr_matrix]:
    # Adapted from https://github.com/scverse/scanpy/blob/29e454429bcb41f3150c2516d82a5c4938f124e1/src/scanpy/tools/_utils.py#L17

    if use_rep is None and n_pcs == 0:
        use_rep = "X"
    if use_rep is None:
        if adata.n_vars > sc.settings.N_PCS:
            if "X_pca" in adata.obsm:
                if n_pcs is not None and n_pcs > adata.obsm["X_pca"].shape[1]:
                    raise ValueError(
                        "`X_pca` does not have enough PCs. Rerun "
                        "`sc.pp.pca` with adjusted `n_comps`."
                    )
                X = adata.obsm["X_pca"][:, :n_pcs]
            else:
                raise ValueError(
                    "You need to run `sc.pp.pca` first. Alternatively, "
                    "specify `use_rep`."
                )
        else:
            X = adata.X
    else:
        if use_rep in adata.obsm and n_pcs is not None:
            if n_pcs > adata.obsm[use_rep].shape[1]:
                raise ValueError(
                    f"{use_rep} does not have enough dimensions. Provide"
                    " a representation with equal or more dimensions "
                    "than `n_pcs` or lower `n_pcs`."
                )
            X = adata.obsm[use_rep][:, :n_pcs]
        elif use_rep in adata.obsm and n_pcs is None:
            X = adata.obsm[use_rep]
        elif use_rep == "X":
            X = adata.X
        else:
            msg = (
                f"Did not find {use_rep} in `.obsm.keys()`. "
                "You need to compute it first."
            )
            raise ValueError(msg)
    return X
