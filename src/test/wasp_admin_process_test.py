import json
import os
import pathlib
import shutil
import tempfile

import numpy as np
from obspy import read
from typer.testing import CliRunner

from .testutils import (
    END_TO_END_DIR,
    RESULTS_DIR,
    get_strong_motion_json,
    get_tele_waves_json,
    update_manager_file_locations,
)

runner = CliRunner()


def test_process_tele():
    from wasp.wasp_admin.process import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        tele_waves = get_tele_waves_json()
        new_tele_waves = update_manager_file_locations(
            tele_waves, tempdir, replace_dir=str(RESULTS_DIR / "data")
        )
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir / "20003k7a_cmt_CMT"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "sampling_filter.json",
            tempdir / "sampling_filter.json",
        )
        os.mkdir(pathlib.Path(tempdir) / "P")
        for o, n in zip(tele_waves, new_tele_waves):
            shutil.copyfile(o["file"], n["file"])
        with open(pathlib.Path(tempdir) / "tele_waves.json", "w") as f:
            json.dump(new_tele_waves, f)

        # test base process tele
        result = runner.invoke(
            app,
            [
                "process-all",
                str(tempdir),
                "body",
                "-g",
                str(tempdir / "20003k7a_cmt_CMT"),
            ],
        )
        assert result.exit_code == 0

        np.testing.assert_array_almost_equal(
            read(new_tele_waves[0]["file"])[0].data[:10],
            [
                -38.251144,
                -38.221546,
                -38.199234,
                -38.189747,
                -38.194893,
                -38.213146,
                -38.240723,
                -38.27213,
                -38.30074,
                -38.320812,
            ],
            5,
        )
        np.testing.assert_array_almost_equal(
            read(new_tele_waves[1]["file"])[0].data[:10],
            [
                -5.4753532,
                -5.5141807,
                -5.5615206,
                -5.61239,
                -5.6537414,
                -5.675031,
                -5.6726656,
                -5.6468496,
                -5.599341,
                -5.537252,
            ],
            5,
        )
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_process_missing_file():
    from wasp.wasp_admin.process import app

    result = runner.invoke(
        app,
        [
            "process-all",
            ".",
            "body",
            "-g",
            str(END_TO_END_DIR / "info" / "20003k7a_cmt_CMT"),
        ],
    )
    assert result.exit_code == 1
    assert isinstance(result.exception, FileNotFoundError)
    assert "does not exist!" in str(result.exception)


def test_remove_baseline():
    from wasp.wasp_admin.process import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        strong_motion = get_strong_motion_json()
        new_strong_motion = update_manager_file_locations(
            strong_motion, tempdir, replace_dir=str(RESULTS_DIR)
        )
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir / "20003k7a_cmt_CMT"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "sampling_filter.json",
            tempdir / "sampling_filter.json",
        )
        os.mkdir(pathlib.Path(tempdir) / "data")
        os.mkdir(pathlib.Path(tempdir) / "data" / "STR")
        for o, n in zip(strong_motion, new_strong_motion):
            shutil.copyfile(o["file"], n["file"])
        with open(pathlib.Path(tempdir) / "strong_motion_waves.json", "w") as f:
            json.dump(new_strong_motion, f)
        # test remove baseline
        result = runner.invoke(
            app,
            [
                "remove-baseline",
                str(tempdir / "data"),
            ],
        )
        assert result.exit_code == 0
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)


def test_shift_match():
    from wasp.wasp_admin.process import app

    tempdir = pathlib.Path(tempfile.mkdtemp())
    try:
        tele_waves = get_tele_waves_json(all=True)
        new_tele_waves = update_manager_file_locations(
            tele_waves, tempdir, replace_dir=str(RESULTS_DIR)
        )
        shutil.copyfile(
            END_TO_END_DIR / "info" / "20003k7a_cmt_CMT", tempdir / "20003k7a_cmt_CMT"
        )
        shutil.copyfile(
            RESULTS_DIR / "NP1" / "sampling_filter.json",
            tempdir / "sampling_filter.json",
        )
        os.mkdir(pathlib.Path(tempdir) / "data")
        os.mkdir(pathlib.Path(tempdir) / "data" / "P")
        os.mkdir(pathlib.Path(tempdir) / "data" / "SH")
        for o, n in zip(tele_waves, new_tele_waves):
            shutil.copyfile(o["file"], n["file"])
        with open(pathlib.Path(tempdir) / "tele_waves.json", "w") as f:
            json.dump(new_tele_waves, f)

        # test base process tele
        result = runner.invoke(
            app,
            [
                "shift-match",
                str(tempdir),
                "body",
                str(tempdir / "20003k7a_cmt_CMT"),
            ],
        )
        result = runner.invoke(
            app,
            [
                "shift-match",
                str(tempdir),
                "body",
                str(tempdir / "20003k7a_cmt_CMT"),
                "-o",
                "manual",
            ],
        )

        assert result.exit_code == 0
    finally:
        print("Cleaning up test directory.")
        shutil.rmtree(tempdir)
