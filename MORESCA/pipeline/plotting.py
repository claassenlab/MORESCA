import logging
from pathlib import Path
from typing import Optional, Union

import gin
import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData

from MORESCA.utils import store_config_params

log = logging.getLogger(__name__)


@gin.configurable(denylist=["sample_id"])
def plotting(
    adata: AnnData,
    apply: bool,
    umap: bool = True,
    path: Union[Path, str] = Path("figures"),
    inplace: bool = True,
    sample_id: Optional[str] = None,
) -> Optional[AnnData]:
    """
    Create plots for an AnnData object.

    Args:
        adata: An AnnData object containing the gene expression data.
        apply: Whether to create plots or not.
        umap: Whether to plot the UMAP or not.
        path: The path to the output directory for the plots.
        inplace: Whether to perform the differential gene expression analysis in-place or return a modified copy of the AnnData object.
        sample_id: Sample ID for subfolder creation in the output directory. Not gin-configurable.

    Returns:
        If `inplace` is True, returns None. Otherwise, returns a modified copy of the AnnData object.
    """
    # TODO: Check before merging if we changed adata
    if not inplace:
        adata = adata.copy()

    store_config_params(
        adata=adata,
        analysis_step=plotting.__name__,
        apply=apply,
        params={
            key: val for key, val in locals().items() if key not in ["adata", "inplace"]
        },
    )

    if not apply:
        return None

    log.info("Creating plots.")

    path = Path(path)

    # Make subfolder if a sample ID is passed (analysis of multiple samples)
    if sample_id:
        path = path / f"{sample_id}/"
    path.mkdir(parents=True, exist_ok=True)

    if umap:
        log.debug("Creating UMAP plot.")
        sc.pl.umap(adata=adata, show=False)
        plt.savefig(Path(path, "umap.png"))
        plt.close()

    if not inplace:
        return adata
