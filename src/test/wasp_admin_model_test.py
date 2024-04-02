from copy import deepcopy
import json
import os
import pathlib
import pytest
import shutil
import tempfile
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
                "run",
                str(tempdir),
                "auto_model",
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
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "segments_data.json", tempdir / "segments_data.json"
        )
        result = runner.invoke(
            app,
            [
                "run",
                str(tempdir),
                "auto_model",
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


@pytest.mark.skipif(
    os.getenv("CI_REGISTRY") is not None,
    reason="Build runner does not have the resources to run",
)
def test_run_multiple():
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
                "run-multiple",
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
                "run-multiple",
                str(tempdir),
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
        shutil.copytree(RESULTS_DIR / "NP1", tempdir / "data" / "NP1")
        # need to update the paths in the surf_waves.json
        with open(tempdir / "data" / "NP1" / "surf_waves.json") as f:
            surfs = json.load(f)
        updated_surfs = []
        for surf_dict in surfs:
            new_surf = deepcopy(surf_dict)
            new_surf["file"] = surf_dict["file"].replace(
                "LONG", str((tempdir / "LONG").resolve())
            )
        updated_surfs += [new_surf]
        with open(tempdir / "data" / "NP1" / "surf_waves.json", "w") as fw:
            json.dump(updated_surfs, fw)
        shutil.copytree(RESULTS_DIR / "data" / "LONG", tempdir / "LONG")
        result = runner.invoke(
            app,
            [
                "run-multiple",
                str(tempdir / "data" / "NP1"),
                "-c",
                str(tempdir / "config.ini"),
                "-t",
                "surf",
                "-g",
                str(tempdir / "20003k7a_cmt_CMT"),
                "-s",
                "9",
                "-s",
                "13",
            ],
        )
        assert result.exit_code == 0
        assert (tempdir / "data" / "NP3.0").exists()
        with open(tempdir / "data" / "NP3.0" / "segments_data.json") as np3:
            np3_data = json.load(np3)
        assert np3_data["segments"][0]["strike"] == 9.0
        assert (tempdir / "data" / "NP3.1").exists()
        with open(tempdir / "data" / "NP3.1" / "segments_data.json") as np31:
            np31_data = json.load(np31)
        assert np31_data["segments"][0]["strike"] == 13.0
    finally:
        print("Cleaning up test directory.")

        shutil.rmtree(tempdir)
