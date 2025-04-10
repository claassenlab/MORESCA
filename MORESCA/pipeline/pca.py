import inspect
import logging
from typing import Optional, Union

import gin
import numpy as np
from anndata import AnnData
from sklearn.decomposition import PCA

from MORESCA.utils import store_config_params

log = logging.getLogger(__name__)


@gin.configurable
def pca(
    adata: AnnData,
    apply: bool,
    n_comps: Union[int, float] = 50,
    use_highly_variable: bool = True,
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Perform principal component analysis (PCA) on the gene expression data in an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to apply the PCA or not.
        n_comps: The number of principal components to compute. A float is interpreted as the proportion of the total variance to retain.
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

    log.info("Performing PCA.")

    if not inplace:
        adata = adata.copy()

    if use_highly_variable:
        log.debug("Using highly variable genes for PCA.")
        X_data = adata[:, adata.var.highly_variable.values].X
    else:
        log.debug("Using all genes for PCA.")
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
