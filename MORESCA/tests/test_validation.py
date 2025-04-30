import gin
import pytest

from MORESCA.pipeline import *  # noqa
from MORESCA.validation import validate_config


def test_config_validation():
    gin.parse_config_file("test-config.gin")
    validate_config(gin.config_str())


def test_faulty_config(caplog):
    gin.parse_config_file("test-config-faulty.gin")
    with pytest.raises(ValueError):
        validate_config(gin.config_str())
    expected_log = (
        "# Errors for clustering:\n"
        "  - Parameter resolution:\n"
        "    - None or more than one of the following rules validate:\n"
        "      - must be of ['integer', 'float'] type\n"
        "      - must be a list  with items validating following rules\n"
        "        - must be of ['integer', 'float'] type\n"
        "        - min value is 0\n"
        "      - must be of string type\n\n"
        "# Errors for diff_gene_exp:\n"
        "  - Parameter method:\n"
        "    - unallowed value Bombardiro Crocodilo\n\n"
        "# Errors for pca:\n"
        "  - Parameter n_comps:\n"
        "    - None or more than one of the following rules validate:\n"
        "      - must be of integer type\n"
        "      - max value is 1.0\n\n"
        "# Errors for quality_control:\n"
        "  - Parameter hb_threshold:\n"
        "    - None or more than one of the following rules validate:\n"
        "      - max value is 100\n"
        "      - must be of ['boolean', 'string'] type\n\n"
        "  - Parameter min_counts:\n"
        "    - None or more than one of the following rules validate:\n"
        "      - min value is 0\n"
        "      - must be of boolean type"
    )
    assert expected_log in caplog.text
