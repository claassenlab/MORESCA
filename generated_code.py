 # Automatically generated code at 12/03/2023 23:00:31
import argparse
import doubletdetection
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scipy.stats as ss
import sys
import warnings
import yaml

from anndata import AnnData
from pathlib import Path
from yaml.loader import SafeLoader

def is_outlier(adata: AnnData, metric: str, nmads: int) -> pd.Series(dtype=bool):
    M = adata.obs[metric]
    MAD = ss.median_abs_deviation(M)
    outlier = (M < np.median(M) - nmads * MAD) | (np.median(M) + nmads * MAD < M)
    return outlier

adata = sc.read(h5adPath)

# Quality control - calculate QC covariates
adata.obs["n_counts"] = adata.X.sum(1)
adata.obs["log_counts"] = np.log(adata.obs["n_counts"])
adata.obs["n_genes"] = (adata.X > 0).sum(1)

adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.contains(("^RP[SL]"))
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt", "ribo", "hb"],
    percent_top=[20],
    log1p=True,
    inplace=True,
)





sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
gene_stack_lst = []



gene_stack_lst.append(np.zeros_like(a=adata.var_names))
remove = np.stack(gene_stack_lst).sum(axis=0).astype(bool)
keep = np.invert(remove)
adata = adata[:, keep]
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=None)
sc.pp.log1p(adata)
sc.pp.normalize_total(adata, target_sum=None)
adata.raw = adata
sc.pp.highly_variable_genes(adata, flavor='seurat_v3', 
            n_top_genes=3000, layer='counts')
sc.pp.scale(adata)
sc.pp.pca(adata, n_comps=50, use_highly_variable=True)

