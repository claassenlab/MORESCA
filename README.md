



[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-3109/)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![codecov](https://codecov.io/gh/claassenlab/MORESCA/branch/main/graph/badge.svg?token=WHUCNFSPJF)](https://codecov.io/gh/claassenlab/MORESCA)

# MORESCA (MOdular and REproducible Single-Cell Analysis)

This repository provides a template  on standardized scRNA-seq analysis using Python and the Scanpy library. All parameters of the workflow are controlled with single config file.

## Usage

### Setting up the environment

Clone the repository 

    git clone git@github.com:claassenlab/MORESCA.git
    
Change into the directory

    cd MORESCA

Create a virtual environment using Conda with Python version >=3.10

    conda create -n <envName> python=3.10

Activate the environment

    conda activate <envName>

Install MORESCA using:

    pip install -e .

This creates a symbolic link, making changes to the code basis instantanious.

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
    doublet_removal = False
    outlier_removal = False
    min_genes = 200
    min_cells = 10
    n_genes_by_counts = None
    mt_threshold = 50
    rb_threshold = 10
    hb_threshold = 2
    remove_mt = False
    remove_rb = False
    remove_hb = False
    remove_custom_genes = None
normalization:
    method = "PFlog1pPF"
feature_selection:
    method = "seurat"
    number_features = 2000
scaling:
    apply = True
    max_value = None
pca:
    apply = True
    n_comps = 50
    use_highly_variable = True
batch_effect_correction:
    method = "harmony"
    batch_key = None
neighborhood_graph:
    n_neighbors = 15
    n_pcs = None
clustering:
    method = "leiden"
    resolution = 1.0
diff_gene_exp:
    method = "wilcoxon"
    groupby = "leiden_r1.0"
    use_raw = True
    tables = False
  ```
  
The following values of the parameters are currently possible

| Parameter | Values 
| - | -
| **quality_control** 
| doublet_removal | *bool* |
| doublet_removal | *bool* |
| min_genes | *int*, *null* | 
| min_cells| *int*, *null* |
| mt_threshold| *float*, *null* |
| rb_threshold| *float*, *null* |
| hb_threshold| *float*, *null* |
| remove_mt| *float*, *null* |
| remove_rb| *float*, *null* |
| remove_hb| *float*, *null* |
| remove_custom_genes| *list(str)*, *null* |
| **normalization**
| method| *log1pCP10k*, *log1PF*, *PFlog1pPF*, *pearson_residuals*, *null*|
| **feature_selection**
| method| *seurat*, *seurat_v3*, *pearson_residuals*, *anti_correlation*, *null*|
| number_features| *int*, *null* |
| **scaling**
| apply| *bool* |
| max_value| *int*, *float* |
| **batch_effect_correction**
| method| *harmony*, *null* |
| batch_key| Not implemented / *null* |
| **neighborhood_graph**
| n_neighbors| *int* |
| n_pcs| *int*, *null* |
| **clustering**
| method| *str*, *null* |
| resolution| *float*|
| **diff_gene_exp**
| method| *wilcoxon*, *logreg*, *t-test*, *t-test_overestim_var* |
| groupy| *str* |
| use_raw| *bool* |
| tables| *bool* |

### Code generator

After deciding for a suitable pipeline and specific parameters, you can create a Python file which reflects the exact step in a minimal fashion. 


## Contributing

For contribution purposes, you should install MORESCA in dev mode:

    pip install -e .[dev]

This additionally install `flake8`, `Black` and `pylint`, which we use for formatting and code style control. Please run these before you commit new code.
Note: This will be made mandatory by using pre-commit hooks.
