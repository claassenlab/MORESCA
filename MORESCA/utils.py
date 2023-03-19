from pathlib import Path
from typing import Optional, Union

from anndata import AnnData


def remove_cells_by_pct_counts(
    adata: AnnData,
    genes: str,
    threshold: Optional[Union[int, float, str, bool]],
    inplace: bool = True,
    save: bool = False,
) -> Optional[AnnData]:
    if not inplace:
        adata = adata.copy()

    if genes not in ["mt", "rb", "hb"]:
        ValueError(
            f"{genes} is not selectable. Accepted values are ['mt', 'rb'', 'hb']"
        )

    match threshold:
        case threshold if isinstance(threshold, (int, float)) and not isinstance(
            threshold, bool
        ):
            if genes == "rb":
                adata = adata[adata.obs[f"pct_counts_{genes}"] > threshold, :]
            else:
                adata = adata[adata.obs[f"pct_counts_{genes}"] < threshold, :]
        case "auto":
            # Todo: Implement automatic selection of threshold.
            ValueError(f"Auto selection for {genes}_threshold not implemented.")
        case False | None:
            print("No mitochondrial filter applied.")
        case _:
            ValueError("Error.")

    if save:
        if isinstance(save, Path | str):
            adata.write(save)
        else:
            # Todo: Revisit this.
            pass

    if not inplace:
        return adata


# Todo: Is this the best way to do it? Manipulating the list inplace feels liek a gotcha.
def remove_genes(gene_lst: list, rmv_lst: list, gene_key: Optional[bool]) -> None:
    match gene_key:
        case True:
            rmv_lst.append(gene_lst)
        case False | None:
            pass
        case _:
            ValueError(f"Invalid choice for remove_{gene_key}.")
