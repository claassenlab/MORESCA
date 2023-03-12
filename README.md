
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-3109/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# scRNAseq-Template

This repository provides a template  on standardized scRNA-seq analysis using Python and the Scanpy library. All parameters of the workflow are controlled with single YAML file.

## Usage

Create a virtual environment using Conda with Python version 3.10:

`conda create -n <envName> python=3.10`

Install the needed packages using the requirements.txt file:

`pip install -r requirements.txt`


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


| Flag | Type | Description | Default |
| - | -  | - | - |
| -d, --data | Path | Path to the h5ad file | *data/adata_raw.h5ad*
| -p, --parameters | Path | Path to the YAML file | *parameters.yml* |
| -v, --verbose | Boolean | If set, prints to output | *False* |
| -f, --figures | Boolean | If set, figures will be generated | *False* |

The following example executes the template with the h5ad file example_data.h5ad, the parameter file example_param.yml and both prints and figures.

```python template.py -d example_data.h5ad -p example_param.yml -v -f```


## Todo

- Better figures (matplotlib.pyplot.subplot_mosaic looks interesting)
- Make all needed methods aware of batch_keys
- Auto-modes where possible
- ..?