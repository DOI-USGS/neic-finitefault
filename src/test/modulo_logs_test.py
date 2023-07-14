import logging
import pathlib
import shutil
import tempfile

from wasp.modulo_logs import add_console_handler, close_log, create_log


def test_create_log():
    tempdir = tempfile.mkdtemp()
    try:
        logger = create_log("test_log_name", pathlib.Path(tempdir) / "testlog.log")
        assert logger.level == 20
        assert logger.name == "test_log_name"
    finally:
        shutil.rmtree(tempdir)


def test_add_console_handler():
    logger = logging.Logger("test_logger")
    updated_logger = add_console_handler(logger, logging.DEBUG)
    assert logger.handlers[0].level == 10


def test_close_log():
    logger = logging.Logger("test_logger")
    logger.addHandler(logging.StreamHandler())
    close_log(logger)
