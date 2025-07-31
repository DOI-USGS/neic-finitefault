# stdlib imports
import configparser
import pathlib
from typing import Union

PROJECT_DIRECTORY = pathlib.Path(__file__).parent.parent.parent
CONFIG_PATH = PROJECT_DIRECTORY / "config.ini"


def _validate_config(config: configparser.ConfigParser, config_path: pathlib.Path):
    """Validate the config file

    :param config: The config object
    :type config: configparser.ConfigParser
    :param config_path: The path to the config file
    :type config_path: pathlib.Path
    :raises ValueError: If a required section and/or option is missing
    """
    errors = []
    options = [
        "code_path",
        "surf_gf_bank",
        "modelling",
        "get_near_gf",
        "compute_near_gf",
        "info",
        "cartopy_files",
    ]
    if "PATHS" not in config.sections():
        errors += [
            (
                "'PATHS' must be a section of the config file "
                f"and include the options:\n{options}"
            )
        ]
    else:
        for option in options:
            if option not in config.options("PATHS"):
                errors += [f"Option in PATHS section '{option}' must be included"]
    if len(errors):
        error_str = "\n".join(errors)
        raise ValueError(f"Error parsing config at {config_path}:\n{error_str}")


def read_config(config_path: pathlib.Path = CONFIG_PATH) -> configparser.ConfigParser:
    """Parse and validate the config file

    :param config_path: The path to the config file, defaults to CONFIG_PATH
    :type config_path: pathlib.Path, optional
    :return: The config parser object
    :rtype: configparser.ConfigParser
    """
    config = configparser.ConfigParser()
    config.read(config_path)
    _validate_config(config=config, config_path=config_path)
    return config


if __name__ == "__main__":
    config = read_config()
    print(config.sections())
    # this will include an empty DEFAULT section,
    # which usefulness is described here:
    # https://www.enfoldsystems.com/software/proxy/docs/4.0/configuringmanually.html#the-default-section
    for section, options in config.items():
        print("Options for section {section}:")
        for key, value in options.items():
            print(key, value)  #'\t{key} = {value}')
    # grabbing values
    print(config["PATHS"]["code_path"])
