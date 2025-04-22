import logging
import warnings
from typing import Optional

import gin
import hotspot
import scanpy as sc
import triku as tk
from anndata import AnnData
from scipy.sparse import csc_matrix

from MORESCA.utils import store_config_params

try:
    from anticor_features.anticor_features import get_anti_cor_genes

    anti_cor_import_error = False
except ImportError:
    anti_cor_import_error = True
    warnings.warn(
        "Could not import anticor_features,\
        install it using 'pip install anticor-features'"
    )

log = logging.getLogger(__name__)


@gin.configurable
def feature_selection(
    adata: AnnData,
    apply: bool,
    method: Optional[str] = "seurat",
    species: Optional[str] = "hsapiens",
    number_features: Optional[int] = None,
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Perform feature selection on an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to apply the feature selection steps or not.
        method: The feature selection method to use. Available options are:
            - "seurat": Use Seurat's highly variable genes method.
            - "seurat_v3": Use Seurat v3's highly variable genes method.
            - "analytical_pearson": Use analytical Pearson residuals for feature selection.
            - "anti_correlation": Use anti-correlation method for feature selection (currently only implemented for human data).
            - "triku": Use Triku to select variable genes.
            - "hotspot": Use Hotspot to select variable genes.
        species: The species for which to perform the anti-correlation feature selection..
        number_features: The number of top features to select (only applicable for certain methods).
        inplace: Whether to perform the feature selection steps in-place or return a modified copy of the AnnData object.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.

    Raises:
        ValueError: If an invalid feature selection method is provided.
    """

    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=feature_selection.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    log.info("Performing feature selection.")

    match method:
        case "seurat":
            log.debug("Using Seurat's highly variable genes method.")
            sc.pp.highly_variable_genes(adata, flavor=method)
        case "seurat_v3":
            log.debug("Using Seurat v3's highly variable genes method.")
            sc.pp.highly_variable_genes(
                adata, flavor=method, n_top_genes=number_features, layer="counts"
            )
        case "analytical_pearson":
            log.debug("Using analytical Pearson residuals for feature selection.")
            sc.experimental.pp.highly_variable_genes(
                adata,
                flavor="pearson_residuals",
                n_top_genes=number_features,
                layer="counts",
            )
        case "anti_correlation":
            log.debug("Using anti-correlation method for feature selection.")

            anti_cor_table = get_anti_cor_genes(
                adata.X.T,
                adata.var.index.tolist(),
                species=species,
                pre_remove_pathways=[],
            )
            anti_cor_table.fillna(value=False, axis=None, inplace=True)
            adata.var["highly_variable"] = anti_cor_table.selected.copy()
        case "triku":
            # Implementation as shown in
            # # https://github.com/theislab/atlas-feature-selection-benchmark/blob/b89fc0f66747062e6e1b4b35bd392b27ad035295/bin/method-triku.py#L15
            log.debug("Using Triku for feature selection.")
            adata_copy = adata.copy()
            adata_copy.X = adata_copy.layers["counts"].copy()

            # Avoid warning about adata.X already being log-transformed
            try:
                del adata_copy.uns
            except AttributeError:
                pass

            sc.pp.filter_cells(adata_copy, min_genes=50)
            sc.pp.filter_genes(adata_copy, min_cells=10)
            sc.pp.normalize_total(adata_copy)

            sc.pp.log1p(adata_copy)
            sc.pp.pca(adata_copy)
            sc.pp.neighbors(
                adata_copy,
                metric="cosine",
                n_neighbors=int(0.5 * len(adata_copy) ** 0.5),
            )
            tk.tl.triku(adata_copy, n_features=number_features)

            adata.var["highly_variable"] = adata_copy.var["highly_variable"]
            del adata_copy

        case "hotspot":
            # Implementation as shown in
            # https://github.com/theislab/atlas-feature-selection-benchmark/blob/b89fc0f66747062e6e1b4b35bd392b27ad035295/bin/method-hotspot.py#L16

            log.debug("Using Hotspot for feature selection.")
            adata_copy = adata.copy()
            adata_copy.X = adata_copy.layers["counts"].copy().astype(int)

            # Avoid warning about adata.X already being log-transformed
            try:
                del adata_copy.uns
            except AttributeError:
                pass

            adata_copy.obs["total_counts"] = adata_copy.X.sum(axis=1)
            adata_copy.layers["counts"] = csc_matrix(adata_copy.X.copy())

            sc.pp.normalize_total(adata_copy)
            sc.pp.log1p(adata_copy)
            sc.pp.scale(adata_copy)
            sc.tl.pca(adata_copy)

            hs = hotspot.Hotspot(
                adata_copy,
                layer_key="counts",
                model="danb",
                latent_obsm_key="X_pca",
                umi_counts_obs_key="total_counts",
            )
            hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
            hs_results = hs.compute_autocorrelations()

            selected_features = (
                hs_results.loc[hs_results["FDR"] < 0.05]
                .sort_values("Z", ascending=False)
                .head(number_features)
                .index
            )

            adata.var["highly_variable"] = adata_copy.var.index.isin(selected_features)
            del adata_copy

        case False | None:
            # TODO: Should this be a warning?
            log.debug("No feature selection applied.")
            return None
        case _:
            raise ValueError(
                f"Selected feature selection method {method} not available."
            )

    if not inplace:
        return adata
