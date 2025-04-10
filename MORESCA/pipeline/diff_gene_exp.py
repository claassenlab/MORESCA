import logging
import warnings
from pathlib import Path
from typing import Literal, Optional, Union

import gin
import pandas as pd
import scanpy as sc
from anndata import AnnData

from MORESCA.utils import store_config_params

log = logging.getLogger(__name__)


@gin.configurable(denylist=["sample_id"])
def diff_gene_exp(
    adata: AnnData,
    apply: bool,
    method: str = "wilcoxon",
    groupby: str = "leiden_r1.0",
    use_raw: Optional[bool] = False,
    layer: Optional[str] = "counts",
    corr_method: Literal["benjamini-hochberg", "bonferroni"] = "benjamini-hochberg",
    tables: Optional[Union[Path, str]] = Path("results/"),
    inplace: bool = True,
    sample_id: Optional[str] = None,
) -> Optional[AnnData]:
    """
    Perform differential gene expression analysis on an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to perform differential gene expression analysis or not.
        method: The differential gene expression analysis method to use. Available options are:
            - "wilcoxon": Use the Wilcoxon rank-sum test.
            - "t-test": Use the t-test.
            - "logreg": Use logistic regression.
            - "t-test_overestim_var": Use the t-test with overestimated variance.
        groupby: The key in `adata.obs` that identifies the groups for comparison.
        use_raw: Whether to use the raw gene expression data or not.
        layer: The layer in `adata.layers` to use for the differential gene expression analysis.
        corr_method: The method to use for multiple testing correction.
        tables: The path to the output directory for the differential expression tables.
        inplace: Whether to perform the differential gene expression analysis in-place or return a modified copy of the AnnData object.
        sample_id: Sample ID for subfolder creation in the output directory. Not gin-configurable.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.

    Note:
        - The result tables are saved as Excel files if `tables` is True.
        - Only genes with adjusted p-values less than 0.05 are included in the result tables.
    """

    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=diff_gene_exp.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    log.info("Determining differentially expressed genes.")

    with warnings.catch_warnings():
        warnings.simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

        # Todo: Should "logreg" be the default?
        match method:
            case method if method in {
                "wilcoxon",
                "t-test",
                "logreg",
                "t-test_overestim_var",
            }:
                key_added = f"{groupby}_{method}"
                log.debug(
                    f"Using method {method}; storing results in adata.uns['{key_added}']."
                )
                sc.tl.rank_genes_groups(
                    adata=adata,
                    groupby=groupby,
                    method=method,
                    corr_method=corr_method,
                    use_raw=use_raw,
                    key_added=key_added,
                    layer=layer,
                )

                dedf_leiden = sc.get.rank_genes_groups_df(
                    adata=adata, group=None, key=key_added
                )

                dedf_leiden.drop("pvals", axis=1, inplace=True)
                # Todo: Should we keep all genes, e.g., for later visualization?
                dedf_leiden = dedf_leiden[dedf_leiden["pvals_adj"] < 0.05]

                if tables:
                    if isinstance(tables, str):
                        tables = Path(tables)
                    if sample_id:
                        tables = Path(tables) / f"{sample_id}/"
                    log.debug(f"Saving differential expression tables to {tables}.")
                    tables.mkdir(parents=True, exist_ok=True)
                    with pd.ExcelWriter(
                        path=f"{tables}/dge_{key_added}.xlsx"
                    ) as writer:
                        for cluster_id in dedf_leiden.group.unique():
                            df_sub_cl = dedf_leiden[
                                dedf_leiden.group == cluster_id
                            ].copy()
                            df_sub_cl.to_excel(writer, sheet_name=f"c{cluster_id}")

            case False | None:
                log.debug("No DGE performed.")
                return None

    if not inplace:
        return adata
