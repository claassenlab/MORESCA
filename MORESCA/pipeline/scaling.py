import logging
from typing import Optional, Union

import gin
import scanpy as sc
from anndata import AnnData

from MORESCA.utils import store_config_params

log = logging.getLogger(__name__)


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

    log.info("Scaling data.")

    sc.pp.scale(adata, max_value=max_value)

    if not inplace:
        return adata
