import logging
from pathlib import Path

import scanpy as sc
from anndata import AnnData

log = logging.getLogger(__name__)


def load_data(data_path) -> AnnData:
    """
    Load data from a specified file path and return an AnnData object.

    Args:
        data_path: The path to the data file.

    Returns:
        An AnnData object containing the loaded data.

    Raises:
        ValueError: If the file format is unknown.

    Note:
        Currently supports loading of '.h5ad', '.loom', and '.h5' file formats.
    """

    if isinstance(data_path, str):
        data_path = Path(data_path)
    if data_path.is_dir():
        # Todo: Implement this for paths.
        pass
    file_extension = data_path.suffix

    log.info(f"Loading data from {data_path} with file extension {file_extension}.")

    match file_extension:
        case ".h5ad":
            adata = sc.read_h5ad(data_path)
        case ".loom":
            adata = sc.read_loom(data_path)
        case ".h5":
            adata = sc.read_10x_h5(data_path)
        case _:
            log.debug(
                f"Unknown file format: {file_extension}. Trying to read with `sc.read`."
            )
            try:
                adata = sc.read(data_path)
            except ValueError:
                raise ValueError(f"Unknown file format: {file_extension}")
    adata.var_names_make_unique()
    return adata
