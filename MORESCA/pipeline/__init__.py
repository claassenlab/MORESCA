import inspect
import logging
import warnings
from pathlib import Path
from typing import List, Literal, Optional, Tuple, Union

import doubletdetection
import gin
import hotspot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scipy.stats as ss
import triku as tk
from anndata import AnnData
from scipy.sparse import csc_matrix
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

from MORESCA.plotting import plot_qc_vars
from MORESCA.utils import (
    choose_representation,
    remove_cells_by_pct_counts,
    remove_genes,
    store_config_params,
)

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
