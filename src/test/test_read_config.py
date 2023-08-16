import configparser
import pathlib

from wasp.read_config import read_config


def test_read_config():
    # should fail if the config is not found or empty
    try:
        config = read_config()
        raise Exception("This should have failed!")
    except ValueError as e:
        assert str(e).endswith(
            "'PATHS' must be a section of the config file and include the options:\n"
            "['code_path', 'surf_gf_bank', 'modelling', 'get_near_gf', 'compute_near_gf', 'info', 'cartopy_files']"
        )

    # should path if the config file is as expected
    config = read_config(
        config_path=pathlib.Path(__file__).parent / "data" / "config.ini"
    )
    assert isinstance(config, configparser.ConfigParser)
