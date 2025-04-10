import logging
from typing import Optional

import gin
import scanpy as sc
from anndata import AnnData

from MORESCA.utils import store_config_params

log = logging.getLogger(__name__)


@gin.configurable
def umap(adata: AnnData, apply: bool, inplace: bool = True) -> Optional[AnnData]:
    """
    Run UMAP on an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to run UMAP or not.
        inplace: Whether to run the UMAP in-place or return a modified copy of the AnnData object.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.
    """
    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=umap.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    log.info("Running UMAP.")

    sc.tl.umap(adata=adata)

    if not inplace:
        return adata
