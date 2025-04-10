import logging
from typing import Optional

import gin
import scanpy as sc
from anndata import AnnData

from MORESCA.utils import store_config_params

log = logging.getLogger(__name__)


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

    log.info("Computing neighborhood graph.")

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
