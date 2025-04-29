import logging
from typing import Optional, Union

from cerberus import Validator

from MORESCA.validation.config_schema import CONFIG_SCHEMA

log = logging.getLogger(__name__)


def _convert_string_to_appropriate_type(
    value: str,
) -> Optional[Union[int, float, str, bool]]:
    """
    Convert a string to its appropriate type.

    Args:
        value: The string to convert.

    Returns:
        The converted value.
    """
    # Numerical values
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass

    # List and tuple values
    if value.startswith("[") and value.endswith("]"):
        return [
            _convert_string_to_appropriate_type(v.strip())
            for v in value[1:-1].split(",")
        ]
    elif value.startswith("(") and value.endswith(")"):
        return tuple(
            _convert_string_to_appropriate_type(v.strip())
            for v in value[1:-1].split(",")
        )

    # None and boolean values
    if value.capitalize() == "None":
        return None
    elif value.capitalize() == "False":
        return False
    elif value.capitalize() == "True":
        return True
    else:
        return value


def _parse_gin_config_str(config_str: str) -> dict:
    """
    Parse the gin config string into a dictionary.

    Args:
        config_str: Configuration string to parse.

    Returns:
        Dictionary representation of the configuration.
    """
    config_dict = {}
    for line in config_str.split("\n"):
        if line.strip() and not line.startswith("#"):
            key, value = line.split(" = ")
            section, parameter = key.split(".")
            section = section.strip()
            parameter = parameter.strip()
            if section not in config_dict:
                config_dict[section] = {}
            value = _convert_string_to_appropriate_type(
                value.strip().strip('"').strip("'")
            )
            config_dict[section][parameter] = value
    return config_dict


def _pretty_print_message_list(message_list: list, indent: int = 0) -> str:
    """
    Pretty print a list of messages.

    Args:
        message_list: List of messages to print.
        indent: Indentation level for pretty printing.

    Returns:
        Formatted string of the messages.
    """
    if len(message_list) == 0:
        return ""
    else:
        message_str = ""
        i = 0
        while i < len(message_list):
            if message_list[i] == "none or more than one rule validate":
                message_str += (
                    f"{' ' * indent}- None or more than one "
                    "of the following rules validate:\n"
                )
                i += 1
                assert isinstance(message_list[i], dict)
                for lst in message_list[i].values():
                    for value in lst:
                        if isinstance(value, dict):
                            message_str += (
                                f"{' ' * (indent + 2)}- must be a list  "
                                "with items validating following rules\n"
                            )
                            for sub_list in value.values():
                                message_str += _pretty_print_message_list(
                                    sub_list, indent=indent + 4
                                )
                        else:
                            message_str += f"{' ' * (indent + 2)}- {value}\n"
            else:
                message_str += f"{' ' * indent}- {message_list[i]}\n"
            i += 1
    return message_str


def _pretty_print_errors(validator) -> str:
    """
    Pretty print the errors from the validator.

    Args:
        validator: The Validator object containing the errors.
    """
    errors_str = "Config validation errors:\n"
    for step, error_list in validator.errors.items():
        errors_str += f"# Errors for {step}:\n"
        for error in error_list:
            for param, message in error.items():
                # TODO: Make message prettier
                errors_str += f"  - Parameter {param}:\n{_pretty_print_message_list(message, indent=4)}\n"
    return errors_str


def validate_config(config_str: str) -> None:
    """
    Validate the gin configuration string.

    Args:
        config_str: Configuration string to validate.

    Raises:
        ValueError: If the configuration is invalid.
    """
    config = _parse_gin_config_str(config_str)

    # TODO: Can we work with gin.config._CONFIG dictionary directly?
    # Make it prettier, see https://github.com/google/gin-config/issues/154

    validator = Validator()
    valid_config = validator.validate(config, CONFIG_SCHEMA)
    if not valid_config:
        log.warning("Configuration validation failed.")
        log.error(_pretty_print_errors(validator))
        raise ValueError("Configuration validation failed. See logs for details.")
    else:
        log.info("Configuration validation passed.")
