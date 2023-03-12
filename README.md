[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-3109/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# scRNAseq-Template

This repository provides a template  on standardized scRNA-seq analysis using Python and the Scanpy library. All parameters of the workflow are controlled with single YAML file.

## Usage

### Setting up the environment

Create a virtual environment using Conda with Python version 3.10:

`conda create -n <envName> -f environment.yml`

Install the needed packages using the requirements.txt file:

The overall folder structure looks like this:

```
project
├── data
│   └── data_raw.h5ad
├── figures
│   ├── postQC
│   └── preQC
│      └── highest_expr_genes.pdf
├── parameters.yml
├── results
│   ├── adata_processed.h5ad
│   └── dge_leiden_r1.0_wilcoxon.xlsx
└── template.py
```

By default, ```template.py``` expects the data in ```H5AD``` format to be in ```data```. The two folders ```figures``` and ```results``` are generated on the fly if they don't exist yet.

Currently, the script will perform the most common operations from doublet removal to DEG analysis of found clusters. If you want to apply ambient RNA correction beforehand, you need to run this separately.

### Calling the template

| Flag | Type | Description | Default |
| - | -  | - | - |
| -d, --data | Path | Path to the h5ad file | *data/adata_raw.h5ad*
| -p, --parameters | Path | Path to the YAML file | *parameters.yml* |
| -v, --verbose | Boolean | If set, prints to output | *False* |
| -f, --figures | Boolean | If set, figures will be generated | *False* |

The following example executes the template with the h5ad file example_data.h5ad, the parameter file example_param.yml and enables both print-statements and figures.

```python template.py -d example_data.h5ad -p example_param.yml -v -f```


### Using the YAML

By default, the used parameter file looks like this:

``` yaml
# parameters.yml
Info:
  sample_key: null
QC:
  doublet_removal: False
  min_genes: 200
  min_cells: 10
  normalization: PFlog1pPF
  feature_selection: seurat_v3
  number_features: 3000
  scale: True
  n_genes_by_counts: null
  mt_threshold: null
  rb_threshold: null
  hb_threshold: null
  remove_mt: False
  remove_rb: False
  remove_hb: False
  remove_custom_genes : null
BatchEffectCorrection:
  method: null
  batch_key: null
NeighborhoodGraph:
  n_neighbors: 15
  n_pcs: null
Clustering:
  method: leiden
  resolution: 1.0
DiffGeneExp:
  method: wilcoxon
```
The following values of the parameters are currently possible


| Parameter | Values 
| - | -
| **INFO** 
| sample_key | *str* |
| **QC** 
| doublet_removal | *bool* |
| min_genes | *int*, *null* | 
| min_cells| *int*, *null* |
| normalization| *log1pCP10k*, *log1PF*, *PFlog1pPF*, *pearson_residuals*, *null*|
| number_features| *seurat*, *seurat_v3*, *pearson_residuals*, *anti_correlation*, *null* |
| feature_selection| *int*, *null* |
| scale| *bool* |
| mt_threshold| *float*, *null* |
| rb_threshold| *float*, *null* |
| hb_threshold| *float*, *null* |
| remove_mt| *float*, *null* |
| remove_rb| *float*, *null* |
| remove_hb| *float*, *null* |
| remove_custom_genes| *list(str)*, *null* |
| **BatchEffectCorrection**
| method| *harmony*, *null* |
| batch_ley| Not implemented / *null* |
| **NeighborhoodGraph**
| n_neighbors| *int* |
| n_pcs| *int*, *null* |
| **Clustering**
| method| *str*, *null* |
| resolution| *float*|
| **DiffGeneExp**
| method| *wilcoxon*, *logreg*, *t-test*, *t-test_overestim_var* |

## Todo

- Better figures (matplotlib.pyplot.subplot_mosaic looks interesting)
- Make all needed methods aware of batch_keys
- Auto-modes where possible
- ..?
