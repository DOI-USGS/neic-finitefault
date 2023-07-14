# -*- coding: utf-8 -*-
"""Module for logger management
"""


import logging


def create_log(name, log_dir):
    """Create a log file with specified name and location
    
    :param name: name of log file
    :param log_dir: location of log file
    :type name: string
    :type log_dir: string
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_dir)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
            '%(asctime)s - %(filename)s - %(lineno)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


def add_console_handler(log, level=logging.INFO):
    """
    """
    ch = logging.StreamHandler()
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    ch.setLevel(level)
    log.addHandler(ch)
    return log


def close_log(log):
    """Create a log file with specified name and location

    :param name: name of log file
    :type name: string
    """
    for i in list(log.handlers):
        log.removeHandler(i)
        i.flush()
        i.close()
    return
