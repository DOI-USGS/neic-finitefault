import pathlib
import shutil
import tempfile
from glob import glob
from unittest import mock

from typer.testing import CliRunner

from .testutils import END_TO_END_DIR, RESULTS_DIR

runner = CliRunner()


@mock.patch("wasp.data_acquisition.acquisition", return_value=None)
def test_acquire(p1):
    from wasp.wasp_admin.data import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir / "20003k7a_cmt_CMT"
        )
        result = runner.invoke(
            app,
            [str(tempdir), str(tempdir / "20003k7a_cmt_CMT")],
        )
    finally:
        print("Cleaning up test directory.")

        shutil.rmtree(tempdir)


def test_acquire_bad_input():
    from wasp.wasp_admin.data import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT",
            tempdir / "20003k7a_cmt_CMT",
        )
        result = runner.invoke(
            app,
            [
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-d",
                "bad_input",
            ],
        )
        assert result.exit_code == 1
        assert (
            result.stdout
            == "'bad_input' is not in the allowed data type list: ['strong', 'tele'].\n"
        )
    finally:
        print("Cleaning up test directory.")

        shutil.rmtree(tempdir)


def test_fill_dicts():
    from wasp.wasp_admin.data import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT",
            tempdir / "20003k7a_cmt_CMT",
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "sampling_filter.json",
            tempdir / "sampling_filter.json",
        )
        result = runner.invoke(
            app,
            [
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-d",
                "tele_body",
            ],
        )
        assert result.exit_code == 0
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_fill_dicts_bad_input():
    from wasp.wasp_admin.data import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT",
            tempdir / "20003k7a_cmt_CMT",
        )
        result = runner.invoke(
            app,
            [
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-d",
                "bad_input",
            ],
        )
        assert result.exit_code == 1
        assert isinstance(result.exception, ValueError)
        assert (
            str(result.exception)
            == "'bad_input' is not in the allowed data type list: ['cgps', 'gps', 'insar', 'strong_motion', 'surf_tele', 'tele_body']."
        )
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_fill_dicts_missing_file():
    from wasp.wasp_admin.data import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT",
            tempdir / "20003k7a_cmt_CMT",
        )
        result = runner.invoke(
            app,
            [
                str(tempdir),
                str(tempdir / "20003k7a_cmt_CMT"),
                "-d",
                "tele_body",
            ],
        )
        assert result.exit_code == 1
        assert isinstance(result.exception, FileNotFoundError)
        assert "does not exist!" in str(result.exception)
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)
