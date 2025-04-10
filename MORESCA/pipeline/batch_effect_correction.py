import logging
from typing import Optional

import gin
import scanpy.external as sce
from anndata import AnnData

from MORESCA.utils import store_config_params

log = logging.getLogger(__name__)


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

    log.info("Performing batch effect correction.")

    match method:
        case "harmony":
            log.debug("Using Harmony for batch effect correction.")
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
            log.debug("No batch effect correction applied.")
            return None
        case _:
            raise ValueError("Invalid choice for batch effect correction method.")

    if not inplace:
        return adata
