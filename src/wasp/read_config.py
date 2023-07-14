# stdlib imports
import configparser
import pathlib
import os

folder = os.path.dirname(__file__)
folder = os.path.dirname(folder)
CONFIG_PATH = os.path.join(os.path.dirname(folder), 'config.ini')
print(CONFIG_PATH)


def read_config():
    config = configparser.ConfigParser()
    config.read(CONFIG_PATH)
    return config


if __name__ == '__main__':
    config = read_config()
    print(config.sections())
    # this will include an empty DEFAULT section,
    # which usefulness is described here:
    # https://www.enfoldsystems.com/software/proxy/docs/4.0/configuringmanually.html#the-default-section
    for section, options in config.items():
        print('Options for section {section}:')
        for key, value in options.items():
            print(key, value)#'\t{key} = {value}')
    # grabbing values
    print(config['PATHS']['code_path'])
