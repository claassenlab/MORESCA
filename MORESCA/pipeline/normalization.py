import logging
from typing import Optional

import gin
import numpy as np
import scanpy as sc
from anndata import AnnData

from MORESCA.utils import remove_genes, store_config_params

log = logging.getLogger(__name__)


@gin.configurable
def normalization(
    adata: AnnData,
    apply: bool,
    method: Optional[str] = "log1pPF",
    remove_mt: Optional[bool] = False,
    remove_rb: Optional[bool] = False,
    remove_hb: Optional[bool] = False,
    inplace: bool = True,
) -> Optional[AnnData]:
    """
    Normalize gene expression data in an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to apply the normalization steps or not.
        method: The normalization method to use. Available options are:
            - "log1pCP10k": Normalize total counts to 10,000 and apply log1p transformation.
            - "log1pPF": Normalize counts per cell to median of total counts and apply log1p transformation.
            - "PFlog1pPF": Normalize counts per cell to median of total counts, apply log1p transformation, and normalize again using the median of total counts.
            - "analytical_pearson": Normalize using analytical Pearson residuals.
        remove_mt: Whether to remove mitochondrial genes or not.
        remove_rb: Whether to remove ribosomal genes or not.
        remove_hb: Whether to remove hemoglobin genes or not.
        inplace: Whether to perform the normalization steps in-place or return a modified copy of the AnnData object.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.

    Raises:
        ValueError: If an invalid normalization method is provided.
    """

    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=normalization.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    log.info("Normalizing total counts per cell.")

    match method:
        case "log1pCP10k":
            log.debug("Normalizing total counts to 10,000 and applying log1p.")
            sc.pp.normalize_total(adata, target_sum=10e4)
            sc.pp.log1p(adata)
        case "log1pPF":
            log.debug(
                "Normalizing counts per cell to median of total counts and applying log1p."
            )
            sc.pp.normalize_total(adata, target_sum=None)
            sc.pp.log1p(adata)
        case "PFlog1pPF":
            log.debug(
                "Normalizing counts per cell to median of total counts, applying log1p, and normalizing again."
            )
            sc.pp.normalize_total(adata, target_sum=None)
            sc.pp.log1p(adata)
            sc.pp.normalize_total(adata, target_sum=None)
        case "analytical_pearson":
            log.debug("Normalizing using analytical Pearson residuals.")
            sc.experimental.pp.normalize_pearson_residuals(adata)
        case None | False:
            log.debug("No normalization applied.")
            return None
        case _:
            raise ValueError(f"Normalization method {method} not available.")

    mt_genes = adata.var_names.str.contains("(?i)^MT-")
    rb_genes = adata.var_names.str.contains("(?i)^RP[SL]")
    hb_genes = adata.var_names.str.contains("(?i)^HB[^(P)]")

    gene_stack_lst = []

    remove_genes(gene_lst=mt_genes, rmv_lst=gene_stack_lst, gene_key=remove_mt)
    remove_genes(gene_lst=rb_genes, rmv_lst=gene_stack_lst, gene_key=remove_rb)
    remove_genes(gene_lst=hb_genes, rmv_lst=gene_stack_lst, gene_key=remove_hb)

    # Add zero array in case all three selection are not selected.
    gene_stack_lst.append(np.zeros_like(a=adata.var_names))
    remove = np.stack(gene_stack_lst).sum(axis=0).astype(bool)
    keep = np.invert(remove)
    adata = adata[:, keep]

    if not inplace:
        return adata
