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
