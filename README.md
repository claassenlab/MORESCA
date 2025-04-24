[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3119/)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/release/python-3128/)
[![Python 3.13](https://img.shields.io/badge/python-3.13-blue.svg)](https://www.python.org/downloads/release/python-3130/)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-red)](https://github.com/astral-sh/ruff)
[![codecov](https://codecov.io/gh/claassenlab/MORESCA/branch/main/graph/badge.svg?token=WHUCNFSPJF)](https://codecov.io/gh/claassenlab/MORESCA)
[![Python package](https://github.com/claassenlab/MORESCA/actions/workflows/python-package.yml/badge.svg)](https://github.com/claassenlab/MORESCA/actions/workflows/python-package.yml)
[![PyPI - Version](https://img.shields.io/pypi/v/moresca)](https://pypi.org/project/moresca/)
[![readthedocs](https://readthedocs.org/projects/MORESCA/badge/?version=latest)](https://moresca.readthedocs.io)


# MORESCA (MOdular and REproducible Single-Cell Analysis)

This repository provides a template  on standardized scRNA-seq analysis using [Python](https://www.python.org/) and the [Scanpy](https://scanpy.readthedocs.io/) library. All parameters of the workflow are controlled with a single config file.

## Usage

### Installation

We strongly recommend to install MORESCA into a virtual environment. Here, we use Conda:

    conda create -n <env_name> python=3.12
    conda activate <env_name>

Then, simply install MORESCA with pip:

    pip install moresca

> [!IMPORTANT]
> If you want to use Python 3.13 on MacOS, make sure to use GCC>=16. This is required for compiling scikit-misc. See [this discussion](https://stackoverflow.com/questions/48174684/fortran-codes-wont-compile-on-mac-with-gfortran) for advice.

### Calling the template

| Flag | Type | Description | Default |
| - | -  | - | - |
| -d, --data | Path | Path to the h5ad file. | *data/adata_raw.h5ad* |
| -o, --output | Path | Path to the output folder for the processed data. | *results* |
| -p, --parameters | Path | Path to the config file. | *config.gin* |
| -v, --verbose | Flag | Generate verbose logs. | |

By default, ```moresca``` expects the data in ```H5AD``` format to be in ```data```. The ```output``` directory as well as figure folders specified in the config are generated on the fly if they don't exist yet.

Currently, the script will perform the most common operations from doublet removal to DEG analysis of found clusters. If you want to apply ambient RNA correction beforehand, you need to run this separately.

The following example executes the pipeline with the h5ad file ```example_data.h5ad``` and the parameter file ```config.gin```, saving the output in the folder ```results```.

    moresca -d example_data.h5ad -o results -p config.gin

### Using the config.gin

By default, the used parameter file looks like this:

```yaml
# config.gin
quality_control:
    apply = True
    doublet_removal = True
    outlier_removal = True
    min_genes = 200
    min_counts = None
    max_counts = None
    min_cells = 10
    max_genes = None
    mt_threshold = 15
    rb_threshold = 10
    hb_threshold = 1
    figures = "figures/"
    pre_qc_plots = True
    post_qc_plots = True

normalization:
    apply = True
    method = "log1pPF"
    remove_mt = False
    remove_rb = False
    remove_hb = False

feature_selection:
    apply = True
    species = "hsapiens"
    method = "seurat_v3"
    number_features = 2000

scaling:
    apply = True
    max_value = None

pca:
    apply = True
    n_comps = 100
    use_highly_variable = True

batch_effect_correction:
    apply = False
    method = "harmony"
    batch_key = None

neighborhood_graph:
    apply = True
    n_neighbors = 30
    n_pcs = None
    metric = "cosine"

clustering:
    apply = True
    method = "leiden"
    resolution = 1.0

diff_gene_exp:
    apply = True
    method = "wilcoxon"
    groupby = "leiden_r1.0"
    use_raw = False
    layer = "unscaled"
    tables = None

umap:
    apply = True

plotting:
    apply = True
    umap = True
    path = "figures/"
```

The following values of the parameters are currently possible

> [!NOTE]
> The parameter ```species``` is only used for the feature selection based on anti-correlation. For available values check [here](https://biit.cs.ut.ee/gprofiler/page/organism-list). It defaults to ```hsapiens```.

| Parameter | Values | Description |
| - | - | - |
| **quality_control** | | |
| apply | *bool* | Whether to apply the quality control steps or not. |
| doublet_removal | *bool* | Whether to perform doublet removal or not. |
| outlier_removal | *bool* | Whether to remove outliers or not. |
| min_genes | *int*, *float*, *bool*, *None* | The minimum number of genes required for a cell to pass quality control. |
| min_counts | *int*, *float*, *bool*, *None* | The minimum total counts required for a cell to pass quality control. |
| max_counts | *int*, *float*, *bool*, *None* | The maximum total counts allowed for a cell to pass quality control. |
| min_cells | *int*, *float*, *bool*, *None* | The minimum number of cells required for a gene to pass quality control. |
| max_genes | *int*, *float*, *str*, *bool*, *None* | The maximum number of genes allowed for a cell to pass quality control. |
| mt_threshold | *int*, *float*, *str*, *bool*, *None* | The threshold for the percentage of counts in mitochondrial genes (maximum). |
| rb_threshold | *int*, *float*, *str*, *bool*, *None* | The threshold for the percentage of counts in ribosomal genes (minimum). |
| hb_threshold | *int*, *float*, *str*, *bool*, *None* | The threshold for the percentage of counts in hemoglobin genes (maximum). |
| figures | *str*, *Path*, *None* | The path to the output directory for the quality control plots. |
| pre_qc_plots | *bool*, *None* | Whether to generate plots of QC covariates before quality control or not. |
| post_qc_plots | *bool*, *None* | Whether to generate plots of QC covariates after quality control or not. |
| **normalization** | | |
| apply | *bool* | Whether to apply the normalization steps or not. |
| method | *log1pCP10k*, *log1pPF*, *PFlog1pPF*, *analytical_pearson*, *None*, *False* | The normalization method to use. |
| remove_mt | *bool*, *None* | Whether to remove mitochondrial genes or not. |
| remove_rb | *bool*, *None* | Whether to remove ribosomal genes or not. |
| remove_hb | *bool*, *None* | Whether to remove hemoglobin genes or not. |
| **feature_selection** | | |
| apply | *bool* | Whether to apply the feature selection steps or not. |
| method | *seurat*, *seurat_v3*, *analytical_pearson*, *anti_correlation*, *triku*, *hotspot*, *None*, *False* | The feature selection method to use. |
| species | *str* | Species of the data. Only used if feature_selection=anti_correlation |
| number_features | *int*, *None* | The number of top features to select (only applicable for certain methods). |
| **scaling** | | |
| apply | *bool* | Whether to apply the scaling step or not. |
| max_value | *int*, *float*, *None* | The maximum value to which the data will be scaled. If None, the data will be scaled to unit variance. |
| **pca** | | |
| apply | *bool* | Whether to apply the PCA or not. |
| n_comps | *int*, *float* | The number of principal components to compute. A float is interpreted as the proportion of the total variance to retain. |
| use_highly_variable | *bool* | Whether to use highly variable genes for PCA computation. |
| **batch_effect_correction** | | |
| apply | *bool* | Whether to apply the batch effect correction or not. |
| method | *harmony*, *None*, *False* | The batch effect correction method to use. |
| batch_key | *str* | The key in `adata.obs` that identifies the batches. |
| **neighborhood_graph** | | |
| apply | *bool* | Whether to compute the neighborhood graph or not. |
| n_neighbors | *int* | The number of neighbors to consider for each cell. |
| n_pcs | *int*, *None* | The number of principal components to use for the computation. |
| metric | *str* | The distance metric to use for computing the neighborhood graph. |
| **clustering** | | |
| apply | *bool* | Whether to perform clustering or not. |
| method | *leiden*, *phenograph*, *None*, *False* | The clustering method to use. |
| resolution | *float*, *int*, *list*, *tuple*, *auto* | The resolution parameter for the clustering method. Can be a single value or a list of values. |
| **diff_gene_exp** | | |
| apply | *bool* | Whether to perform differential gene expression analysis or not. |
| method | *wilcoxon*, *t-test*, *logreg*, *t-test_overestim_var* | The differential gene expression analysis method to use. |
| groupby | *str* | The key in `adata.obs` that identifies the groups for comparison. |
| use_raw | *bool*, *None* | Whether to use the raw gene expression data or not. |
| layer | *str*, *None* | The layer in `adata.layers` to use for the differential gene expression analysis. |
| corr_method | *benjamini-hochberg*, *bonferroni* | The method to use for multiple testing correction. |
| tables | *str*, *Path*, *None* | The path to the output directory for the differential expression tables. |
| **umap** | | |
| apply | *bool* | Whether to run UMAP or not. |
| **plotting** | | |
| apply | *bool* | Whether to create plots or not. |
| umap | *bool* | Whether to plot the UMAP or not. |
| path | *str*, *Path* | The path to the output directory for the plots. |

## Contributing

For contribution purposes, you should clone MORESCA from GitHub and install it in dev mode:

    git clone git@github.com:claassenlab/MORESCA.git
    cd MORESCA
    pip install -e ".[dev]"
    pre-commit install

This additionally installs `ruff` and `pytest`, which we use for formatting and code style control. Please run these before you commit new code.
Additionally, it will set up a pre-commit hook to run `ruff`.
