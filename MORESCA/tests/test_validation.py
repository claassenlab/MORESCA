import gin
import pytest

from MORESCA.pipeline import *  # noqa
from MORESCA.validation import validate_config


def test_config_validation():
    gin.parse_config_file("test-config.gin")
    validate_config(gin.config_str())


def test_faulty_config():
    gin.parse_config_file("test-config-faulty.gin")
    with pytest.raises(ValueError):
        validate_config(gin.config_str())
