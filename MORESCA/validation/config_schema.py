_QC_POSITIVE_NUMBER_PARAMETERS = {
    param: {
        "oneof": [
            {"type": ["integer", "float"], "min": 0},
            {"type": "boolean", "allowed": [False]},
        ],
        "nullable": True,
    }
    for param in ["min_genes", "min_counts", "max_counts", "min_cells", "max_genes"]
}
_QC_PERCENTAGE_PARAMETERS = {
    param: {
        "oneof": [
            {"type": ["integer", "float"], "min": 0, "max": 100},
            {"type": ["boolean", "string"], "allowed": [False, "auto"]},
        ],
        "nullable": True,
    }
    for param in ["mt_threshold", "rb_threshold", "hb_threshold"]
}

_QC_SCHEMA = dict(
    {
        "apply": {"type": "boolean", "required": True},
        "doublet_removal": {"type": "boolean"},
        "outlier_removal": {"type": "boolean"},
        "figures": {"type": "string", "nullable": True},
        "pre_qc_plots": {"type": "boolean", "nullable": True},
        "post_qc_plots": {"type": "boolean", "nullable": True},
    },
    **_QC_POSITIVE_NUMBER_PARAMETERS,
    **_QC_PERCENTAGE_PARAMETERS,
)

_NORMALIZATION_SCHEMA = dict(
    {
        "apply": {"type": "boolean", "required": True},
        "method": {
            "type": ["string", "boolean"],
            "nullable": True,
            "allowed": [
                "log1pCP10k",
                "log1pPF",
                "PFlog1pPF",
                "analytical_pearson",
                False,
            ],
        },
    },
    **{
        param: {"type": "boolean", "nullable": True}
        for param in ["remove_mt", "remove_rb", "remove_hb"]
    },
)

_FEATURE_SELECTION_SCHEMA = {
    "apply": {"type": "boolean", "required": True},
    "method": {
        "type": ["string", "boolean"],
        "nullable": True,
        "allowed": [
            "seurat",
            "seurat_v3",
            "analytical_pearson",
            "anti_correlation",
            "triku",
            "hotspot",
            False,
        ],
    },
    "species": {
        "type": "string",
        "nullable": True,
    },  # TODO: make required for anti_correlation method
    "number_features": {"type": "integer", "nullable": True, "min": 1},
}

_CLUSTERING_SCHEMA = {
    "apply": {"type": "boolean", "required": True},
    "method": {
        "type": ["string", "boolean"],
        "nullable": True,
        "allowed": ["leiden", "phenograph", False],
    },
    "resolution": {
        "oneof": [
            {"type": ["integer", "float"], "min": 0},
            {"type": "list", "schema": {"type": ["integer", "float"], "min": 0}},
            {"type": "string", "allowed": ["auto"]},
        ]
    },
}

_DIFF_GENE_EXP_SCHEMA = {
    "apply": {"type": "boolean", "required": True},
    "method": {
        "type": "string",
        "allowed": ["wilcoxon", "t_test", "logreg", "t-test_overestim_var"],
    },
    "groupby": {"type": "string"},
    "use_raw": {"type": "boolean", "nullable": True},
    "layer": {"type": "string", "nullable": True},
    "corr_method": {"type": "string", "allowed": ["benjamini-hochberg", "bonferroni"]},
    "tables": {"type": "string", "nullable": True},
}

CONFIG_SCHEMA = {
    "quality_control": {"type": "dict", "schema": _QC_SCHEMA},
    "normalization": {"type": "dict", "schema": _NORMALIZATION_SCHEMA},
    "feature_selection": {"type": "dict", "schema": _FEATURE_SELECTION_SCHEMA},
    "scaling": {
        "type": "dict",
        "schema": {
            "apply": {"type": "boolean", "required": True},
            "max_value": {"type": ["integer", "float"], "nullable": True},
        },
    },
    "pca": {
        "type": "dict",
        "schema": {
            "apply": {"type": "boolean", "required": True},
            "n_comps": {
                "oneof": [
                    {"type": "integer", "min": 1},
                    {"type": "float", "min": 0.0, "max": 1.0},
                ]
            },
            "use_highly_variable": {"type": "boolean"},
        },
    },
    "batch_effect_correction": {
        "type": "dict",
        "schema": {
            "apply": {"type": "boolean", "required": True},
            "method": {
                "type": ["string", "boolean"],
                "nullable": True,
                "allowed": ["harmony", False],
            },
            "batch_key": {
                "type": ["string"],
                "nullable": True,
            },  # TODO: only nullable if apply or method is False
        },
    },
    "neighborhood_graph": {
        "type": "dict",
        "schema": {
            "apply": {"type": "boolean", "required": True},
            "n_neighbors": {"type": "integer", "min": 1},
            "n_pcs": {"type": "integer", "min": 1, "nullable": True},
            "metric": {"type": "string"},
        },
    },
    "clustering": {"type": "dict", "schema": _CLUSTERING_SCHEMA},
    "diff_gene_exp": {"type": "dict", "schema": _DIFF_GENE_EXP_SCHEMA},
    "umap": {
        "type": "dict",
        "schema": {"apply": {"type": "boolean", "required": True}},
    },
    "plotting": {
        "type": "dict",
        "schema": {
            "apply": {"type": "boolean", "required": True},
            "umap": {"type": "boolean"},
            "path": {"type": "string"},
        },
    },
}
