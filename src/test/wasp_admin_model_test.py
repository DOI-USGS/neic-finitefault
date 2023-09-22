import json
import os
import pathlib
import shutil
import tempfile
from glob import glob
from unittest import mock

import numpy as np
from obspy import read
from typer.testing import CliRunner

from .testutils import DATA_DIR, END_TO_END_DIR, HOME, RESULTS_DIR

runner = CliRunner()


@mock.patch(target="wasp.inversion_chen_new.automatic_usgs", return_value=None)
def test_run(p1):
    from wasp.wasp_admin.model import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    with open(HOME / "fortran_code" / "gfs_nm" / "long" / "low.in", "r") as f:
        file_data = f.read()
    with open(HOME / "fortran_code" / "gfs_nm" / "long" / "low.in", "a") as f:
        if "fd_bank" not in file_data:
            f.write(
                f"\n{str((HOME  / 'fortran_code' / 'gfs_nm'/ 'long'/'fd_bank').resolve())}"
            )
    try:
        shutil.copytree(END_TO_END_DIR / "data", tempdir / "data")

        with open(DATA_DIR / "config.ini") as f:
            config = f.read().replace("/home/user/neic-finitefault", str(HOME))
        with open(tempdir / "config.ini", "w") as wf:
            wf.write(config)

        # test no cmt
        result = runner.invoke(
            app,
            [
                str(tempdir),
                "usgs_model",
                "-d",
                str(tempdir / "data"),
                "-c",
                str(tempdir / "config.ini"),
                "-t",
                "body",
            ],
        )
        assert result.exit_code == 1
        assert (
            str(result.exception)
            == "Either gcmt_tensor_file or qcmt_tensor_file must be defined."
        )

        # test file does not exist
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir / "20003k7a_cmt_CMT"
        )
        result = runner.invoke(
            app,
            [
                str(tempdir),
                "usgs_model",
                "-d",
                str(tempdir / "data"),
                "-c",
                str(tempdir / "config.ini"),
                "-t",
                "body",
                "-g",
                str(tempdir / "20003k7a_cmt_CMT"),
            ],
        )
        assert result.exit_code == 1
        assert "does not exist!" in str(result.exception)

        # test run
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "segments_data.json", tempdir / "segments_data.json"
        )
        result = runner.invoke(
            app,
            [
                str(tempdir),
                "usgs_model",
                "-d",
                str(tempdir / "data"),
                "-c",
                str(tempdir / "config.ini"),
                "-t",
                "body",
                "-g",
                str(tempdir / "20003k7a_cmt_CMT"),
            ],
        )
        assert result.exit_code == 0
        assert (tempdir / "20150916225432" / "ffm.0" / "NP1").exists()
        assert (tempdir / "20150916225432" / "ffm.0" / "data").exists()
        assert (tempdir / "20150916225432" / "ffm.0" / "logs").exists()
    finally:
        print("Cleaning up test directory.")

        shutil.rmtree(tempdir)
