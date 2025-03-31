[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3119/)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/release/python-3128/)
[![Python 3.13](https://img.shields.io/badge/python-3.13-blue.svg)](https://www.python.org/downloads/release/python-3130/)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-red)](https://github.com/astral-sh/ruff)
[![codecov](https://codecov.io/gh/claassenlab/MORESCA/branch/main/graph/badge.svg?token=WHUCNFSPJF)](https://codecov.io/gh/claassenlab/MORESCA)
[![Python package](https://github.com/claassenlab/MORESCA/actions/workflows/python-package.yml/badge.svg)](https://github.com/claassenlab/MORESCA/actions/workflows/python-package.yml)


# MORESCA (MOdular and REproducible Single-Cell Analysis)

This repository provides a template  on standardized scRNA-seq analysis using Python and the Scanpy library. All parameters of the workflow are controlled with single config file.

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
| -d, --data | Path | Path to the h5ad file | *data/adata_raw.h5ad*
| -p, --parameters | Path | Path to the config file | *config.gin* |
| -v, --verbose | Boolean | If set, prints to output | *False* |
| -f, --figures | Boolean | If set, figures will be generated | *False* |

By default, ```template.py``` expects the data in ```H5AD``` format to be in ```data```. The two folders ```figures``` and ```results``` are generated on the fly if they don't exist yet.

Currently, the script will perform the most common operations from doublet removal to DEG analysis of found clusters. If you want to apply ambient RNA correction beforehand, you need to run this separately.

The following example executes the template with the h5ad file example_data.h5ad, the parameter file config.gin and enables both print-statements and figures.

```python template.py -d example_data.h5ad -p config.gin -v -f```

### Using the config.gin

By default, the used parameter file looks like this:

``` yml
# config.gin
quality_control:
    apply = True
    doublet_removal = True
    outlier_removal = True
    min_genes = 200
    min_counts = None
    max_counts = None
    min_cells = 10
    n_genes_by_counts = None
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

| Parameter | Values
| - | -
| **quality_control**
| apply | *bool* |
| doublet_removal | *bool* |
| doublet_removal | *bool* |
| min_genes | *int*, *null* |
| min_cells| *int*, *null* |
| mt_threshold| *float*, *null* |
| rb_threshold| *float*, *null* |
| hb_threshold| *float*, *null* |
| figures| *str*|
| pre_qc_plots | *bool* |
| post_qc_plots | *bool* |
| **normalization**
| method| *log1pCP10k*, *log1PF*, *PFlog1pPF*, *pearson_residuals*, *null*|
| remove_mt| *bool*, *null* |
| remove_rb| *bool*, *null* |
| remove_hb| *bool*, *null* |
| remove_custom_genes| Not implemented |
| **feature_selection**
| apply| *bool* |
| method| *seurat*, *seurat_v3*, *pearson_residuals*, *anti_correlation*, *null*|
| number_features| *int*, *null* |
| **scaling**
| apply| *bool* |
| max_value| *int*, *float* |
| **pca**
| apply| *bool* |
| n_comps| *int*, *float* |
| use_highly_variable| *bool*|
| **batch_effect_correction**
| apply| *bool* |
| method| *harmony*, *null* |
| batch_key| Not implemented / *null* |
| **neighborhood_graph**
| apply| *bool* |
| n_neighbors| *int* |
| n_pcs| *int*, *null* |
| metric| *str* |
| **clustering**
| apply| *bool* |
| method| *str*, *null* |
| resolution| *float*|
| **diff_gene_exp**
| apply| *bool* |
| method| *wilcoxon*, *logreg*, *t-test*, *t-test_overestim_var* |
| groupby| *str* |
| use_raw| *bool* |
| layer| *str*, *null* |
| tables| *bool* |

## Contributing

For contribution purposes, you should install MORESCA in dev mode:

    pip install -e ".[dev]"

This additionally installs `ruff` and `pytest`, which we use for formatting and code style control. Please run these before you commit new code.
Note: This will be made mandatory by using pre-commit hooks.
