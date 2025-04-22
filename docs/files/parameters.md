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
| method | *leiden*, *None*, *False* | The clustering method to use. |
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
