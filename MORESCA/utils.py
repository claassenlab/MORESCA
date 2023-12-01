from pathlib import Path
from typing import Optional, Union
import numpy as np
import scipy.stats as ss
import scanpy as sc
from anndata import AnnData


def is_passing_upper(data, nmads, upper_limit=0):
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
    upper_bound = max(med + nmads * mad, upper_limit)
    return data <= upper_bound


def is_passing_lower(data, nmads, lower_limit):
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
    save: bool = False,
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

    match threshold:
        case threshold if isinstance(threshold, (int, float)) and not isinstance(
            threshold, bool
        ):
            if genes == "rb":
                adata = adata[adata.obs[f"pct_counts_{genes}"] > threshold, :]
            else:
                adata = adata[adata.obs[f"pct_counts_{genes}"] < threshold, :]
        case "auto":
            raise ValueError(f"Auto selection for {genes}_threshold not implemented.")
        case False | None:
            print(f"No {genes} filter applied.")
        case _:
            raise ValueError("Error.")

    if save and isinstance(save, Path | str):
        adata.write(save)
    if not inplace:
        return adata


# Todo: Is this the best way to do it? Manipulating the list inplace feels like a gotcha.
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


def ddqc(adata: AnnData, inplace: bool = True, save: bool = False) -> Optional[AnnData]:
    """
    Perform Data-Driven Quality Control (DDQC) on an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        inplace: Whether to perform the DDQC in-place or return a modified copy of the AnnData object.
        save: Whether to save the modified AnnData object to a file. If True, the AnnData object will be saved.

    Returns:
        If `inplace` is True and `save` is False, returns None. Otherwise, returns a modified copy of the AnnData object.

    Note:
        - The resulting AnnData object will only contain cells that pass the quality control checks.
    """

    if not inplace:
        adata = adata.copy()

    adata_raw = adata.copy()

    adata_raw = adata.copy()
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    adata = adata[adata.obs.pct_counts_mt <= 80, :]
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=2000, layer="counts"
    )
    sc.pp.scale(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50, metric="euclidean")
    sc.tl.leiden(adata, resolution=1.0)

    # Create a numpy array to store the passing information
    n_obs = adata.n_obs
    passed = np.ones(n_obs, dtype=bool)

    # Map clusters to contiguous integer labels
    cluster_idx_map = {
        cluster: np.where(adata.obs["leiden"].values == cluster)[0]
        for cluster in adata.obs["leiden"].unique()
    }

    # Iterate through each cluster
    for indices in cluster_idx_map.values():
        # Subset relevant cluster data
        pct_counts_mt_cluster = adata.obs["pct_counts_mt"].values[indices]
        total_counts_cluster = adata.obs["total_counts"].values[indices]
        n_genes_cluster = adata.obs["n_genes"].values[indices]

        # Apply the upper bound check for pct_counts_mt
        passing_mask_mt = is_passing_upper(pct_counts_mt_cluster, nmads=3)

        # Apply the lower bound check for n_counts and n_genes
        passing_mask_counts = is_passing_lower(
            total_counts_cluster, nmads=3, lower_limit=0
        )
        passing_mask_genes = is_passing_lower(n_genes_cluster, nmads=3, lower_limit=200)

        # Combine the masks
        full_passing_mask = passing_mask_mt & passing_mask_counts & passing_mask_genes

        # Update the 'passed' numpy array using the cluster indices
        passed[indices] = full_passing_mask

    # Subselect 'adata' based on the 'passed' mask
    adata = adata_raw[passed].copy()

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
            key: (False if key == "apply" else None) for key, val in params.items()
        }
    adata.uns[uns_key][analysis_step] = {
        key: (list(val) if isinstance(val, tuple) else val)
        for key, val in params.items()
    }
