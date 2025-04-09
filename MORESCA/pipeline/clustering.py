import logging
from typing import List, Literal, Optional, Tuple, Union

import gin
import numpy as np
import scanpy as sc
from anndata import AnnData
from sklearn.metrics import silhouette_score

from MORESCA.utils import choose_representation, store_config_params

log = logging.getLogger(__name__)


@gin.configurable
def clustering(
    adata: AnnData,
    apply: bool,
    method: str = "leiden",
    resolution: Union[
        float, int, List[Union[float, int]], Tuple[Union[float, int]], Literal["auto"]
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

    log.info("Performing clustering.")

    match method:
        case "leiden":
            log.debug("Using Leiden algorithm for clustering.")
            if (
                not isinstance(resolution, (float, int, list, tuple))
                and resolution != "auto"
            ):
                raise ValueError(f"Invalid type for resolution: {type(resolution)}.")

            if isinstance(resolution, (float, int)):
                log.debug(f"Using single resolution {resolution} for clustering.")
                resolutions = [resolution]
            elif resolution == "auto":
                log.debug("Using auto resolution for clustering.")
                resolutions = [0.25] + list(np.linspace(0.5, 1.5, 11)) + [2.0]
                log.debug(f"Tested resolutions: {[float(r) for r in resolutions]}.")
            else:
                log.debug(f"Using multiple resolutions {resolution} for clustering.")
                resolutions = resolution

            for res in resolutions:
                sc.tl.leiden(
                    adata=adata,
                    resolution=res,
                    flavor="igraph",
                    n_iterations=2,
                    key_added=f"leiden_r{res}",
                    random_state=0,
                )
        case False | None:
            log.debug("No clustering applied.")
            return None
        case _:
            raise ValueError(f"Clustering method {method} not available.")

    # Choose best resolution according to silhouette score
    if len(resolutions) > 1:
        neighbors_params = adata.uns["neighbors"]["params"]
        metric = neighbors_params["metric"]
        use_rep = (
            None if "use_rep" not in neighbors_params else neighbors_params["use_rep"]
        )
        n_pcs = None if "n_pcs" not in neighbors_params else neighbors_params["n_pcs"]

        # Use the representation used for neighborhood graph computation
        X = choose_representation(adata, use_rep=use_rep, n_pcs=n_pcs)

        scores = np.zeros(len(resolutions))

        for i, res in enumerate(resolutions):
            scores[i] = silhouette_score(
                X, labels=adata.obs[f"leiden_r{res}"], metric=metric
            )

        best_res = resolutions[np.argmax(scores)]
        log.debug(f"Best resolution: {best_res}.")
        adata.obs["leiden"] = adata.obs[f"leiden_r{best_res}"]

        adata.uns["MORESCA"]["clustering"]["best_resolution"] = best_res
        adata.uns["MORESCA"]["clustering"]["resolutions"] = resolutions
        adata.uns["MORESCA"]["clustering"]["silhouette_scores"] = scores

    if not inplace:
        return adata
