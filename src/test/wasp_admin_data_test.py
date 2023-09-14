import pathlib
import shutil
import tempfile
from glob import glob
from unittest import mock

from typer.testing import CliRunner

from .testutils import END_TO_END_DIR

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
