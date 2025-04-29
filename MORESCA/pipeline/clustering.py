import logging
from typing import List, Literal, Optional, Tuple, Union

import gin
import numpy as np
import scanpy as sc
import scanpy.external as sce
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
            - "phenograph": Use the Phenograph algorithm for clustering.
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

    method_config = {
        "leiden": (
            sc.tl.leiden,
            {"n_iterations": -1, "flavor": "igraph", "random_state": 0},
        ),
        "phenograph": (
            sce.tl.phenograph,
            {"clustering_algo": "leiden", "n_iterations": -1, "seed": 0},
        ),
    }

    func_to_call, params_to_use = method_config[method]

    match method:
        case "leiden" | "phenograph":
            log.debug(f"Using {method} algorithm for clustering.")
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
                if method == "leiden":
                    params_to_use["resolution"] = res
                    params_to_use["key_added"] = f"{method}_r{res}"
                elif method == "phenograph":
                    params_to_use["resolution_parameter"] = res
                func_to_call(adata, **params_to_use)

                if method == "phenograph":
                    adata.obs.rename(
                        columns={"pheno_leiden": f"phenograph_r{res}"}, inplace=True
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
                X, labels=adata.obs[f"{method}_r{res}"], metric=metric
            )

        best_res = resolutions[np.argmax(scores)]
        log.debug(f"Best resolution: {best_res}.")
        adata.obs[method] = adata.obs[f"{method}_r{best_res}"]

        adata.uns["MORESCA"]["clustering"]["best_resolution"] = best_res
        adata.uns["MORESCA"]["clustering"]["resolutions"] = resolutions
        adata.uns["MORESCA"]["clustering"]["silhouette_scores"] = scores

    if not inplace:
        return adata
