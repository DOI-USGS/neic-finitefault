# -*- coding: utf-8 -*-
"""Module for logger management
"""


import logging
import pathlib
from typing import Union


def create_log(name: str, log_file: Union[pathlib.Path, str]) -> logging.Logger:
    """Create a log file with specified name and location

    :param name: The log name
    :type name: str
    :param log_file: The path to the log file
    :type log_file: Union[pathlib.Path, str]
    :return: _description_
    :rtype: logging.Logger
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_file)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s - %(filename)s - %(lineno)s - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


def add_console_handler(
    log: logging.Logger, level: int = logging.INFO
) -> logging.Logger:
    """Add logging to the console

    :param log: The logger
    :type log: logging.Logger
    :param level: The log level, defaults to logging.INFO
    :type level: int, optional
    :return: The updated logger
    :rtype: logging.Logger
    """
    ch = logging.StreamHandler()
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    ch.setLevel(level)
    log.addHandler(ch)
    return log


def close_log(log: logging.Logger):
    """Create a log file with specified name and location

    :param log: The logger
    :type log: logging.Logger
    """
    for i in list(log.handlers):
        log.removeHandler(i)
        i.flush()
        i.close()
    return
