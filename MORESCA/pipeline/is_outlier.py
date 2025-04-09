import logging

import numpy as np
import pandas as pd
import scipy.stats as ss
from anndata import AnnData

log = logging.getLogger(__name__)


def is_outlier(adata: AnnData, metric: str, nmads: int) -> pd.Series:
    """
    Check if each value in a given metric column of an AnnData object is an outlier.

    Args:
        adata: An AnnData object containing the data.
        metric: The name of the metric column to check for outliers.
        nmads: The number of median absolute deviations (MADs) away from the median to consider a value as an outlier.

    Returns:
        A pandas Series of boolean values indicating whether each value is an outlier or not.
    """
    data = adata.obs[metric]
    med_abs_dev = ss.median_abs_deviation(data)
    return (data < np.median(data) - nmads * med_abs_dev) | (
        np.median(data) + nmads * med_abs_dev < data
    )
